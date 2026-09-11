module BioformatsLoader
using JavaCall
import JavaCall.JNI
using LightXML
using ImageMetadata
using AxisArrays
using URIs
using ImageCore
using Downloads
using Logging

export
    set_id!,
    OMEXMLReader,
    openbytes,
    open_reinterpret,
    set_series!,
    getpixeltype,
    open_img,
    bf_import,
    open_stack,
    metadata,
    arraydata,
    axisnames,
    get_num_sets,
    get_RGB_channel_count,
    get_bf_path,
    default_bf_path,
    is_available,
    download_bioformats!

include("metadata.jl")
include("omexmlreader.jl")

const scheme2importer = Dict{String, Any}()

function add_importer(importer, schemes)
    for s in schemes
        scheme2importer[s] = importer
    end
end

function open_img(oxr::OMEXMLReader, index::Int)
    sx = get_size(oxr, "X")
    sy = get_size(oxr, "Y")
    arr = openbytes(oxr, index)
    return interpret_blob!(oxr, arr, (sx, sy))
end

for F in (:open_stack, :bf_import, :metadata)
    @eval begin
        function $F(filename::String, redirect::IO)
            redirect_stdout(() -> $F(filename), redirect)
        end
        function $F(filename::String, show_output::Bool)
            if show_output
                $F(filename)
            else
                STDOLD = stdout
                local ret
                redirect_stdout()
                try
                    ret = $F(filename)
                finally
                    redirect_stdout(STDOLD)
                end
                ret
            end
        end
    end
end

function open_stack(filename::String)
    OMEXMLReader(open_stack, filename)
end

function interpret_blob!(oxr, arr, sizes)
    arr = reinterpret(getpixeltype(oxr), arr)
    if !islittleendian(oxr)
        @. arr = ntoh(arr)
    end
    return reshape(arr, sizes)
end

function no_colon(sz, subidx)
    ntuple((d)->typeof(subidx[d]) == Colon ? (1:sz[d]) : 
            ((typeof(subidx[d]) <: Number) ? (subidx[d]:subidx[d]) : subidx[d]), length(sz))
end

function sub_indices(sz, subidx)
    CartesianIndices(no_colon(sz, subidx))
end

function clipped_size(sz, subidx)
    return size(sub_indices(sz, subidx))
end


function open_stack(oxr::OMEXMLReader; subidx=nothing, order="CYXZT")
    image_order = get_dimension_order(oxr)
    if isinterleaved(oxr)
        image_order = replace(image_order, "XYC" => "CXY")
    end

    sizes = [get_size(oxr, d) for d in image_order]
    size_dict = Dict(zip(image_order, sizes))

    for d in image_order
        if !in(d, order)
            @assert (size_dict[d] == 1) "Non-singleton dimension \"$d\" is not present in the requested image order \"$order\""
        end
    end

    num_imgs = get_image_count(oxr)
    trail_size=1;
    trail_dim=length(image_order)
    for d=length(image_order):-1:1
        trail_size *= size_dict[image_order[d]]
        if trail_size == num_imgs
            trail_dim = d
            break
        elseif trail_size > num_imgs
            error("size inconsistency")
        end
    end
    raw_size = Tuple(size_dict[d] for d in image_order if in(d, order))

    subidx = let
        if isnothing(subidx)
            [:,:,:,:,:] # 0:(num_imgs - 1)
        else
            sub_dict = Dict(zip(order, subidx))
            Tuple(sub_dict[d] for d in image_order if in(d, order))
        end
    end
    fsubidx = subidx[1:trail_dim-1]
    tsubidx = subidx[trail_dim:end]

    all_idx = sub_indices(raw_size[trail_dim:end], tsubidx)
    LI = LinearIndices(raw_size[trail_dim:end])
    arr = open_reinterpret(oxr, LI[all_idx[1]]-1, fsubidx, raw_size[1:trail_dim-1])
    for ci in all_idx[2:end]
        idx = LI[ci]
        arr_nd = open_reinterpret(oxr, idx-1, fsubidx, raw_size[1:trail_dim-1])
        append!(arr, arr_nd)
    end
    new_raw_size = clipped_size(raw_size, subidx)
    im = interpret_blob!(oxr, arr, (new_raw_size...,))
    im = permutedims(im, [findfirst(c, image_order) for c in order])

    return AxisArray(im, (Symbol(d) for d in order)...)
end

"""
    bf_import(uri::AbstractString; order::AbstractString="TZYXC", squeeze::Bool=false, subset=nothing, subidx=nothing)

Import all images and metadata in the file using Bio-Formats.

The `order` keyword argument is the output order of dimensions, allowing you to
reshuffle the dimensions of the images or skip singleton dimensions. Providing
"CYX" for `order` gives you three-dimensional images in a "channels-first" data
format but will fail with an error if any of the images in the file has a size
greater than one in the `T` or `Z` dimension.

By setting the `squeeze` keyword argument to `true` all singleton dimensions in
the images will be dropped (even if they are in the `order` argument).

The keyword argument `subset` allows to only import one or multiple specific sets (for datasets with multiple sets)
and the keyword argument `subidx` allows to import specific ranges of the data (in the `order` as given by the user).
"""
function bf_import(uri; kwargs...)
    u = URI(uri)
    if u.scheme in keys(scheme2importer)
        scheme2importer[u.scheme](u; kwargs...)
    else
        if length(u.scheme) > 1
            @warn "Unrecognized scheme \"$(u.scheme)\", attempting to open as file"
        end
        bf_import_file(uri; kwargs...)
    end
end

"""
    get_num_sets(filename::AbstractString)

returns the number of series present in the file as given by `filename`. Currently URLs are not allowed to avoid multiple downloads. 
"""
function get_num_sets(filename::AbstractString)
    OMEXMLReader(filename) do oxr
        xml = get_xml(oxr)
        metalst = get_elements_by_tagname(root(xml), "Image")
        return length(metalst)
    end
end

bf_import_file(uri::URI; kwargs...) = bf_import_file(uri.path; kwargs...)
function bf_import_file(filename::AbstractString; subset=nothing, subidx=nothing, order="CYXZT", squeeze=false, gray=false)
    OMEXMLReader(filename) do oxr
        xml = get_xml(oxr)
        images = Array{ImageMeta,1}()

        metalst = get_elements_by_tagname(root(xml), "Image")
        subset = isnothing(subset) ? (1:length(metalst)) : subset
        for i in subset
            set_series!(oxr, i-1)
            img = open_stack(oxr; subidx=subidx, order = order)
            if squeeze
                img = dropdims(img; dims=((i for i in 1:ndims(img) if size(img, i) == 1)...,))
            end
            if gray
                img = as_gray(img)
            end
            properties = xml_to_dict(metalst[i])
            properties[:ImportOrder] = axisnames(img)
            push!(images, ImageMeta(img, properties))
        end

        return images
    end
end

add_importer(bf_import_file, ["file"])

"""
    bf_import_http(url::AbstractString, [ filename::AbstractString ]; kwargs...)

Download and import an image using Bio-Formats.
"""
bf_import_http(url::AbstractString, args...; kwargs...) = bf_import_http(URI(url), args...; kwargs...)

function bf_import_http(url::URI; kwargs...)
    filename = basename(url.path)
    @assert (filename != "") "Could not determine filename from URL"
    bf_import_http(url, filename; kwargs...)
end

function bf_import_http(url::URI, filename::AbstractString; kwargs...)
    mktempdir() do path
        imgpath = joinpath(path, filename)
        Downloads.download(string(url), imgpath)
        bf_import_file(imgpath; kwargs...)
    end
end

add_importer(bf_import_http, ["http", "https"])

function as_gray(arr::A) where A <: AbstractArray{T} where T <: Unsigned
    map(x -> Gray(reinterpret(Normed{T, sizeof(T)*8}, x)), arr)
end
function as_gray(arr::A) where A <: AbstractArray{T} where T <: AbstractFloat
    map(Gray, arr)
end
function as_gray(arr::A) where A <: AbstractArray{Bool}
    map(Gray, arr)
end

function metadata(filename::String)
    OMEXMLReader(filename) do oxr
        xml = get_xml(oxr)
        properties = Dict{String, Any}()
        metalst = get_elements_by_tagname(root(xml), "Image")

        SizeC = get_size(oxr, "C")

        for a = attributes(metalst[1])
            push!(properties, "$(name(a))"=>value(a))
        end
        for a = attributes(find_element(metalst[1], "Pixels"))
            try
                push!(properties, "$(name(a))"=> parse(Float64, value(a)))
            catch ex
                if ex isa ArgumentError
                    push!(properties, "$(name(a))"=>  value(a))
                end
            end
        end

        if length(metalst[1]["Pixels"][1]["Channel"]) == SizeC
            for i in 1:SizeC
                for a in attributes(metalst[1]["Pixels"][1]["Channel"][i])
                    if name(a)== "EmissionWavelength"
                        push!(properties, "channel_$i"=>value(a))
                    end
                end
            end
        end
        properties
    end
end

const levels = ["ALL", "DEBUG", "ERROR", "FATAL", "INFO", "OFF", "TRACE", "WARN"]
function enableLogging(level::String)
    level in levels || error("enableLogging: level argument must be one of "*join(levels,", "))
    dt = @jimport(loci.common.DebugTools)
    Bool(jcall(dt,"enableLogging",jboolean,(JString,),level))
end

"""
    default_bf_path() -> String

Preferred location for `bioformats_package.jar`: `bioformats/bioformats_package.jar`
under the first entry of `Base.DEPOT_PATH`. This is the location
`download_bioformats!` writes to, and (unless overridden, see [`get_bf_path`](@ref))
the location `get_bf_path` resolves to.
"""
default_bf_path() = joinpath(Base.DEPOT_PATH[1], "bioformats", "bioformats_package.jar")

# Legacy location: next to the package itself, filled in by the old
# `Pkg.build`-time download (pre-dating `default_bf_path`/`download_bioformats!`).
_legacy_bf_path() = joinpath(dirname(@__DIR__), "deps", "bioformats_package.jar")

"""
    get_bf_path() -> String

Resolve the path to `bioformats_package.jar`, checked in order:
1. The `BIOFORMATS_JAR_PATH` environment variable, if set and pointing at an
   existing file.
2. [`default_bf_path()`](@ref) (the Julia depot), if the file exists there.
3. The legacy package-local `deps/bioformats_package.jar` path, if it exists
   there (back-compat with the old `Pkg.build`-time download).

Falls back to [`default_bf_path()`](@ref) if none of the above exist, so
callers get a consistent target path to download to. Use [`is_available`](@ref)
to check whether the resolved path actually exists.
"""
function get_bf_path()
    env = get(ENV, "BIOFORMATS_JAR_PATH", "")
    !isempty(env) && isfile(env) && return env
    isfile(default_bf_path()) && return default_bf_path()
    isfile(_legacy_bf_path()) && return _legacy_bf_path()
    return default_bf_path()
end

"""
    is_available() -> Bool

Whether `bioformats_package.jar` currently exists at the path [`get_bf_path`](@ref)
resolves to.
"""
is_available() = isfile(get_bf_path())

"""
    download_bioformats!(; version="6.5.1", dest=default_bf_path(), progress=nothing) -> String

Download `bioformats_package.jar` (GPL v2, published by the Open Microscopy
Environment) from `downloads.openmicroscopy.org` to `dest`, creating parent
directories as needed. Downloads to a temporary `.part` file first and renames
it into place atomically, so an interrupted download never leaves a file
[`is_available`](@ref) would wrongly accept. Returns `dest`.

`progress`, if given, is passed through to `Downloads.download` as the
progress callback (see `Downloads.download`'s `progress` keyword).
"""
function download_bioformats!(; version::AbstractString="6.5.1", dest::AbstractString=default_bf_path(), progress=nothing)
    url = "https://downloads.openmicroscopy.org/bio-formats/$version/artifacts/bioformats_package.jar"
    mkpath(dirname(dest))
    tmp = dest * ".part"
    try
        if progress === nothing
            Downloads.download(url, tmp)
        else
            Downloads.download(url, tmp; progress=progress)
        end
        mv(tmp, dest; force=true)
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
    return dest
end

function set_memory(memory=1024)
    if memory > 0
        # There should be a distinct memory option in JavaCall
        # Find other memory options and delete them
        filter!(opt -> !startswith(opt, "-Xmx"), JavaCall.opts)
        # Add a new option as specified
        JavaCall.addOpts("-Xmx$(memory)M")
    end
    nothing
end

function __init__()
    # Configure static JVM options only. Deliberately does NOT call
    # `get_bf_path()`/`addClassPath` here: `__init__` runs at module/process
    # initialization time, which for an app embedding this package in a
    # PackageCompiler sysimage is *before* the app's own startup code gets a
    # chance to run (e.g. fixing up `Base.DEPOT_PATH` to a location that
    # doesn't depend on how the app was launched — `get_bf_path()`'s default
    # resolves via `Base.DEPOT_PATH[1]`). Resolving the classpath this early
    # can silently freeze in a stale/wrong path. `get_bf_path()` is instead
    # resolved lazily in `init()`, right before the JVM actually launches.
    JavaCall.addOpts("-ea")
    JavaCall.addOpts("-Xrs")
    set_memory()
end

let initialized = Ref(false)
    global init, ensure_init
    function init(;memory::Int=-1, log_level::String="ERROR")
        if !initialized[]
            if !JavaCall.isloaded()
                JavaCall.addClassPath(get_bf_path())
                set_memory(memory)
                JavaCall.init()
                atexit(JavaCall.destroy)
            else
                @warn "JavaCall already initialized"
            end
            initialized[] = true
            enableLogging(log_level) || @warn "Could not enable logging."
        else
            @warn "BioformatsLoader already initialized"
        end
        nothing
    end
    ensure_init() = initialized[] ? nothing : init()
end

end

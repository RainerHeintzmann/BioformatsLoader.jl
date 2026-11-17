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
    series_info,
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
    # Bio-Formats accepts both local file paths and remote URLs (http/https).
    # A local path is not a valid URI in general: parsing it with `URI(...)`
    # throws for any path containing a space or control character ("raw
    # whitespace" — RFC 3986), and even a plain Windows path without spaces
    # risks its drive letter ("C:") being misparsed as a URI scheme. So only
    # attempt URI parsing when `uri` actually looks like `<scheme>://...`;
    # anything else (including every local path) goes straight to
    # `bf_import_file`, exactly as it would have if `URI(uri)` had "succeeded"
    # with an unrecognized/empty scheme.
    if occursin(r"^[a-zA-Z][a-zA-Z0-9+.\-]*://", uri)
        u = URI(uri)
        if u.scheme in keys(scheme2importer)
            return scheme2importer[u.scheme](u; kwargs...)
        end
        @warn "Unrecognized scheme \"$(u.scheme)\", attempting to open as file"
    end
    bf_import_file(uri; kwargs...)
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

"""
    series_info(filename::AbstractString; thumbnails=true, time_budget=20.0) -> Vector{Dict{Symbol,Any}}

Describe every series (OME `Image`) in `filename` without reading its pixel data,
e.g. to let a user choose one from a multi-series file (LIF, CZI, ND2, …). Each
entry has `:index` (1-based, as `bf_import`'s `subset`), `:name` (the OME image
name, may be empty), `:size` (`[X, Y, Z, C, T]`), `:dimension_order` (e.g.
`"XYCZT"`), `:pixel_type` and, with `thumbnails=true`, `:thumbnail`: a
`(width, height, pixels)` named tuple of 8-bit gray values (x fastest) made by
Bio-Formats' `openThumbBytes` from the middle Z plane of the first channel and
time point (or `nothing`). Thumbnails stop once `time_budget` seconds are used
up, so a file with very many (or very large) series still answers quickly.
"""
function series_info(filename::AbstractString; thumbnails::Bool=true, time_budget::Real=20.0)
    OMEXMLReader(filename) do oxr
        xml = get_xml(oxr)
        images = get_elements_by_tagname(root(xml), "Image")
        t0 = time()
        out = Vector{Dict{Symbol,Any}}()
        for (i, img) in enumerate(images)
            set_series!(oxr, i - 1)
            d = Dict{Symbol,Any}(
                :index => i,
                :name => has_attribute(img, "Name") ? attribute(img, "Name") : "",
                :size => [get_size(oxr, c) for c in "XYZCT"],
                :dimension_order => get_dimension_order(oxr),
                :pixel_type => string(getpixeltype(oxr)),
            )
            if thumbnails && time() - t0 < time_budget
                d[:thumbnail] = try
                    _series_thumbnail(oxr)
                catch e
                    @debug "BioformatsLoader: no thumbnail for series $i" exception = e
                    nothing
                end
            end
            push!(out, d)
        end
        return out
    end
end

# 8-bit gray thumbnail of the current series (middle Z plane, first channel and
# time point), contrast-stretched between the 0.5 % and 99.5 % percentiles.
function _series_thumbnail(oxr::OMEXMLReader)
    r = oxr.reader
    z = get_size(oxr, "Z") ÷ 2
    idx = jcall(r, "getIndex", jint, (jint, jint, jint), z, 0, 0)
    w = Int(jcall(r, "getThumbSizeX", jint, ()))
    h = Int(jcall(r, "getThumbSizeY", jint, ()))
    bytes = local_frame() do
        jcall(r, "openThumbBytes", Vector{jbyte}, (jint,), idx)
    end
    T = getpixeltype(oxr)
    vals = reinterpret(T, bytes)
    islittleendian(oxr) || (vals = ntoh.(vals))
    n = w * h
    nrgb = max(1, get_RGB_channel_count(oxr))
    # RGB thumbnails hold `nrgb` samples per pixel: keep the first one
    v = nrgb > 1 && length(vals) >= n * nrgb ?
        (isinterleaved(oxr) ? vals[1:nrgb:n*nrgb] : vals[1:n]) : vals[1:min(n, length(vals))]
    f = Float32.(v)
    length(f) == n || return nothing
    s = sort(filter(isfinite, f))
    isempty(s) && return nothing
    lo = s[max(1, round(Int, 0.005 * length(s)))]
    hi = s[max(1, round(Int, 0.995 * length(s)))]
    scale = hi > lo ? 255f0 / (hi - lo) : 0f0
    px = UInt8[isfinite(x) ? round(UInt8, clamp((x - lo) * scale, 0f0, 255f0)) : 0x00 for x in f]
    return (width=w, height=h, pixels=px)
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
Bio-Formats version that [`download_bioformats!`](@ref) fetches by default (the
newest stable release; its classes still run on Java 8).
"""
const BIOFORMATS_VERSION = "8.5.0"

"""
    bioformats_urls(version=BIOFORMATS_VERSION) -> Vector{String}

Where `bioformats_package.jar` of `version` is downloaded from, tried in this order:
1. OME's Maven repository (`artifacts.openmicroscopy.org`): fixed, versioned file
   names that build tools depend on, so the most stable source;
2. the GitHub release of Bio-Formats (`v<version>`, available since about 6.6);
3. `downloads.openmicroscopy.org` (the link on OME's download page, which has
   been unavailable at times).
"""
bioformats_urls(version::AbstractString=BIOFORMATS_VERSION) = [
    "https://artifacts.openmicroscopy.org/artifactory/ome.releases/ome/bioformats_package/$version/bioformats_package-$version.jar",
    "https://github.com/ome/bioformats/releases/download/v$version/bioformats_package.jar",
    "https://downloads.openmicroscopy.org/bio-formats/$version/artifacts/bioformats_package.jar",
]

# A complete jar (zip) file: starts with a local file header and has the end of
# central directory record near its end. Rejects HTML error pages and truncated downloads.
function _is_complete_jar(path::AbstractString)
    n = filesize(path)
    n > 1_000_000 || return false
    open(path) do io
        read(io, 4) == UInt8[0x50, 0x4b, 0x03, 0x04] || return false
        tail = min(n, 66_000)   # the record (22 bytes) is followed by a comment of at most 64 KB
        seek(io, n - tail)
        buf = read(io, tail)
        eocd = UInt8[0x50, 0x4b, 0x05, 0x06]
        return any(i -> view(buf, i:i+3) == eocd, length(buf)-21:-1:1)
    end
end

"""
    download_bioformats!(; version=BIOFORMATS_VERSION, dest=default_bf_path(), urls=bioformats_urls(version),
                          progress=nothing, log=nothing) -> String

Download `bioformats_package.jar` (GPL v2, published by the Open Microscopy
Environment) to `dest`, creating parent directories as needed. The `urls` are
tried in order (see [`bioformats_urls`](@ref)) until one delivers a complete jar;
if all fail, the error lists every URL with its error. Downloads to a temporary
`.part` file first and renames it into place, so an interrupted download never
leaves a file [`is_available`](@ref) would wrongly accept. Returns `dest`.

`progress`, if given, is passed through to `Downloads.download` as the
progress callback (see `Downloads.download`'s `progress` keyword). `log`, if
given, is called with a message string before each attempt and after a failed one.
"""
function download_bioformats!(; version::AbstractString=BIOFORMATS_VERSION, dest::AbstractString=default_bf_path(),
                              urls=bioformats_urls(version), progress=nothing, log=nothing)
    say(msg) = log === nothing ? @info(msg) : log(msg)
    mkpath(dirname(dest))
    tmp = dest * ".part"
    failures = String[]
    try
        for url in urls
            say("Downloading Bio-Formats $version from $url ...")
            try
                if progress === nothing
                    Downloads.download(url, tmp)
                else
                    Downloads.download(url, tmp; progress=progress)
                end
                _is_complete_jar(tmp) || error("not a complete jar file ($(filesize(tmp)) bytes)")
                mv(tmp, dest; force=true)
                return dest
            catch e
                e isa InterruptException && rethrow()
                msg = sprint(showerror, e)
                push!(failures, "$url: $msg")
                say("  failed: $msg")
                isfile(tmp) && rm(tmp; force=true)
            end
        end
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
    error("Could not download bioformats_package.jar (Bio-Formats $version) from any source:\n  " *
          join(failures, "\n  ") * "\nPlace it by hand at $dest instead (then restart).")
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

# JavaCall (≤ 0.8.1, upstream master as of 2026-11) keeps one JNI environment per
# thread in `JNI.ppenv`, indexed by `Threads.threadid()`, but sizes it with
# `Threads.nthreads()` — the size of the *default* pool only. Since Julia 1.12
# every session also has an interactive thread (`julia -t 1` = 1 default +
# 1 interactive thread, ids 1 and 2), so `threadid()` exceeds `nthreads()` and
# `JNI.attach_threads()` throws a BoundsError right after the JVM was created
# (it looked like "no Java found"). Every later JNI call made on such a thread
# would hit the same BoundsError. Repair: size the table by `maxthreadid()` and
# attach the remaining threads ourselves. No-op where JavaCall got it right.
function _javacall_init()
    try
        JavaCall.init()
    catch e
        # The JVM itself is up (only attaching the other threads failed):
        # repair below. Anything else (no Java, bad options) is a real error.
        JavaCall.isloaded() || rethrow()
        @debug "BioformatsLoader: repairing JavaCall's per-thread JNI table" exception = e
    end
    _attach_all_threads!()
end

function _attach_all_threads!()
    ppenv = JavaCall.JNI.ppenv
    n = Threads.maxthreadid()
    length(ppenv) < n && append!(ppenv, fill(Ptr{JavaCall.JNI.JNIEnv}(C_NULL), n - length(ppenv)))
    attach!() = (tid = Threads.threadid();
                 ppenv[tid] == C_NULL && JavaCall.JNI.attach_current_thread(Ref(ppenv, tid)))
    attach!()   # the calling thread (normally the one that created the JVM)
    # Other threads: best effort. Without JULIA_COPY_STACKS the JVM refuses to
    # attach threads that run Julia tasks on their own stacks ("Unable to attach
    # thread id"), and JULIA_COPY_STACKS=1 can hang Julia 1.12 at start-up. That
    # is fine as long as all JNI calls come from the thread that called `init()`
    # (DeconvGUI/ImageProcessingGUI run them in a JVM worker loop on the main
    # task), so a failure here must not make the JVM unusable.
    failed = Threads.Atomic{Int}(0)
    Threads.@threads :static for _ in 1:Threads.nthreads()   # one per default-pool thread
        try
            attach!()
        catch e
            Threads.atomic_add!(failed, 1)
            @debug "BioformatsLoader: could not attach thread $(Threads.threadid()) to the JVM" exception = e
        end
    end
    failed[] > 0 && @debug "BioformatsLoader: $(failed[]) thread(s) not attached to the JVM; " *
                           "call Bio-Formats only from thread $(Threads.threadid())."
    nothing
end

let initialized = Ref(false)
    global init, ensure_init
    function init(;memory::Int=-1, log_level::String="ERROR")
        if !initialized[]
            if !JavaCall.isloaded()
                JavaCall.addClassPath(get_bf_path())
                set_memory(memory)
                _javacall_init()
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

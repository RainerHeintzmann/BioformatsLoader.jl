xsd_version = "2016-06"
xsd_url = "https://www.openmicroscopy.org/Schemas/OME/$xsd_version/ome.xsd"
@info "Downloading version $xsd_version of ome.xsd from $xsd_url"
download(xsd_url, joinpath(pwd(), "ome.xsd"))

# bioformats_package.jar is GPL-licensed and, unlike the OME-XML schema above,
# is not required for the package to load — it is only needed at runtime, once
# `BioformatsLoader.init()` actually starts the JVM. To avoid silently fetching
# GPL-licensed code on every `Pkg.add`/`Pkg.build`, it is not downloaded here by
# default; set `BIOFORMATS_AUTO_DOWNLOAD=true` to restore the old behavior, or
# call `BioformatsLoader.download_bioformats!()` / set `BIOFORMATS_JAR_PATH` at
# runtime instead.
if get(ENV, "BIOFORMATS_AUTO_DOWNLOAD", "") == "true"
    # (the same version and sources as `BioformatsLoader.download_bioformats!()`)
    version = "8.5.0"
    bfpkg_urls = [
        "https://artifacts.openmicroscopy.org/artifactory/ome.releases/ome/bioformats_package/$version/bioformats_package-$version.jar",
        "https://github.com/ome/bioformats/releases/download/v$version/bioformats_package.jar",
        "https://downloads.openmicroscopy.org/bio-formats/$version/artifacts/bioformats_package.jar",
    ]
    ok = false
    for bfpkg_url in bfpkg_urls
        @info "Downloading version $version of bioformats_package.jar from $bfpkg_url"
        try
            download(bfpkg_url, joinpath(pwd(), "bioformats_package.jar"))
            global ok = true
            break
        catch e
            @warn "Download failed" bfpkg_url exception = e
        end
    end
    ok || error("Could not download bioformats_package.jar from any of: " * join(bfpkg_urls, ", "))
else
    @info "Skipping bioformats_package.jar download (GPL-licensed). Set ENV[\"BIOFORMATS_AUTO_DOWNLOAD\"]=\"true\" before Pkg.build to restore the old behavior, or fetch it at runtime via BioformatsLoader.download_bioformats!() / BIOFORMATS_JAR_PATH."
end

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
    version = "6.5.1"
    bfpkg_url = "https://downloads.openmicroscopy.org/bio-formats/$version/artifacts/bioformats_package.jar"
    @info "Downloading version $version of bioformats_package.jar from $bfpkg_url"
    download(bfpkg_url, joinpath(pwd(), "bioformats_package.jar"))
else
    @info "Skipping bioformats_package.jar download (GPL-licensed). Set ENV[\"BIOFORMATS_AUTO_DOWNLOAD\"]=\"true\" before Pkg.build to restore the old behavior, or fetch it at runtime via BioformatsLoader.download_bioformats!() / BIOFORMATS_JAR_PATH."
end

# Download, unzip, check permissions, and test the bioregion's binary files

This function downloads and unzips the 'bin' folder required to run
certain functions of the `bioregion` package. It also verifies if the
files have the necessary permissions to be executed as programs.
Finally, it tests whether the binary files are running correctly.

## Usage

``` r
install_binaries(
  binpath = "tempdir",
  download_only = FALSE,
  infomap_version = c("2.1.0", "2.6.0", "2.7.1", "2.8.0"),
  verbose = TRUE
)
```

## Arguments

- binpath:

  A `character` string specifying the path to the folder that will host
  the `bin` folder containing the binary files (see Details).

- download_only:

  A `logical` value indicating whether the function should only download
  the `bin.zip` file or perform the entire process (see Details).

- infomap_version:

  A `character` vector or a single `character` string specifying the
  Infomap version(s) to install.

- verbose:

  A `boolean` indicating whether to display progress messages. Set to
  `FALSE` to suppress these messages.

## Value

No return value.

## Details

By default, the binary files are installed in R's temporary directory
(`binpath = "tempdir"`). In this case, the `bin` folder will be
automatically removed at the end of the R session. Alternatively, the
binary files can be installed in the `bioregion` package folder
(`binpath = "pkgfolder"`).

A custom folder path can also be specified. In this case, and only in
this case, `download_only` can be set to `TRUE`, but you must ensure
that the files have the required permissions to be executed as programs.

**In all cases, PLEASE MAKE SURE to update the `binpath` and
`check_install` parameters accordingly in
[netclu_infomap](https://bioRgeo.github.io/bioregion/reference/netclu_infomap.md),
[netclu_louvain](https://bioRgeo.github.io/bioregion/reference/netclu_louvain.md),
and
[netclu_oslom](https://bioRgeo.github.io/bioregion/reference/netclu_oslom.md).**

## Note

Currently, only Infomap versions 2.1.0, 2.6.0, 2.7.1, and 2.8.0 are
available.

## See also

For more details illustrated with a practical example, see the vignette:
<https://biorgeo.github.io/bioregion/articles/a1_install_binary_files.html>.

## Author

Maxime Lenormand (<maxime.lenormand@inrae.fr>)  
Boris Leroy (<leroy.boris@gmail.com>)  
Pierre Denelle (<pierre.denelle@gmail.com>)

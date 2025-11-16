# 1. Installation of the binary files

Some functions (listed below) or at least part of them require binary
files to run.

- [netclu_infomap](https://bioRgeo.github.io/bioregion/reference/netclu_infomap.html)
- [netclu_louvain](https://bioRgeo.github.io/bioregion/reference/netclu_louvain.html)
  (Cpp version)
- [netclu_oslom](https://bioRgeo.github.io/bioregion/reference/netclu_oslom.html)

The function
[install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.html)
downloads and unzips the `bin` folder needed to run these functions. It
also checks if the files have the permissions to be executed as
programs. It finally tests if the binary files are running properly.

## 1. Run install_binaries()

The function
[install_binaries](https://bioRgeo.github.io/bioregion/reference/install_binaries.html)
should be run prior to using
[netclu_infomap](https://bioRgeo.github.io/bioregion/reference/netclu_infomap.html),
the Cpp version of
[netclu_louvain](https://bioRgeo.github.io/bioregion/reference/netclu_louvain.html)
and
[netclu_oslom](https://bioRgeo.github.io/bioregion/reference/netclu_oslom.html)
as follows.

``` r
install_binaries(binpath = "tempdir" ,
                 download_only =  FALSE,
                 infomap_version = c("2.1.0", 
                                     "2.6.0",
                                     "2.7.1",
                                     "2.8.0"))
```

    ## 

    ## 1. Download bin.zip

    ## 

    ## The folder has been successfully downloaded to /tmp/RtmpR44Nyi.

    ## 

    ## 2. Unzip folder

    ## 

    ## The folder has been successfully unzipped in /tmp/RtmpR44Nyi.

    ## 

    ## 3. Check permissions

    ## 

    ## Permission to execute /tmp/RtmpR44Nyi/bin/INFOMAP/2.1.0/infomap_omp_lin as program: denied.

    ## Permission to execute /tmp/RtmpR44Nyi/bin/INFOMAP/2.1.0/infomap_noomp_lin as program: denied.

    ## Permission to execute /tmp/RtmpR44Nyi/bin/INFOMAP/2.6.0/infomap_omp_lin as program: denied.

    ## Permission to execute /tmp/RtmpR44Nyi/bin/INFOMAP/2.6.0/infomap_noomp_lin as program: denied.

    ## Permission to execute /tmp/RtmpR44Nyi/bin/INFOMAP/2.7.1/infomap_omp_lin as program: denied.

    ## Permission to execute /tmp/RtmpR44Nyi/bin/INFOMAP/2.7.1/infomap_noomp_lin as program: denied.

    ## Permission to execute /tmp/RtmpR44Nyi/bin/INFOMAP/2.8.0/infomap_omp_lin as program: denied.

    ## Permission to execute /tmp/RtmpR44Nyi/bin/INFOMAP/2.8.0/infomap_noomp_lin as program: denied.

    ## Permission to execute /tmp/RtmpR44Nyi/bin/LOUVAIN/convert_lin as program: denied.

    ## Permission to execute /tmp/RtmpR44Nyi/bin/LOUVAIN/louvain_lin as program: denied.

    ## Permission to execute /tmp/RtmpR44Nyi/bin/OSLOM/oslom_dir_lin as program: denied.

    ## Permission to execute /tmp/RtmpR44Nyi/bin/OSLOM/oslom_undir_lin as program: denied.

    ## 

    ## Try to change permissions automatically

    ## 

    ## Automatic change of permission succeed, /tmp/RtmpR44Nyi/bin/INFOMAP/2.1.0/infomap_omp_lin can now be executed as progam.

    ## Automatic change of permission succeed, /tmp/RtmpR44Nyi/bin/INFOMAP/2.1.0/infomap_noomp_lin can now be executed as progam.

    ## Automatic change of permission succeed, /tmp/RtmpR44Nyi/bin/INFOMAP/2.6.0/infomap_omp_lin can now be executed as progam.

    ## Automatic change of permission succeed, /tmp/RtmpR44Nyi/bin/INFOMAP/2.6.0/infomap_noomp_lin can now be executed as progam.

    ## Automatic change of permission succeed, /tmp/RtmpR44Nyi/bin/INFOMAP/2.7.1/infomap_omp_lin can now be executed as progam.

    ## Automatic change of permission succeed, /tmp/RtmpR44Nyi/bin/INFOMAP/2.7.1/infomap_noomp_lin can now be executed as progam.

    ## Automatic change of permission succeed, /tmp/RtmpR44Nyi/bin/INFOMAP/2.8.0/infomap_omp_lin can now be executed as progam.

    ## Automatic change of permission succeed, /tmp/RtmpR44Nyi/bin/INFOMAP/2.8.0/infomap_noomp_lin can now be executed as progam.

    ## Automatic change of permission succeed, /tmp/RtmpR44Nyi/bin/LOUVAIN/convert_lin can now be executed as progam.

    ## Automatic change of permission succeed, /tmp/RtmpR44Nyi/bin/LOUVAIN/louvain_lin can now be executed as progam.

    ## Automatic change of permission succeed, /tmp/RtmpR44Nyi/bin/OSLOM/oslom_dir_lin can now be executed as progam.

    ## Automatic change of permission succeed, /tmp/RtmpR44Nyi/bin/OSLOM/oslom_undir_lin can now be executed as progam.

    ## 

    ## 4. Test Infomap

    ## 

    ## Congratulation, you successfully install the 2.1.0 OpenMP version of Infomap!

    ## Congratulation, you successfully install the 2.6.0 OpenMP version of Infomap!

    ## Congratulation, you successfully install the 2.7.1 OpenMP version of Infomap!

    ## Congratulation, you successfully install the 2.8.0 OpenMP version of Infomap!

    ## 

    ## 5. Test Louvain

    ## 

    ## Congratulation, you successfully install the version 0.3 of Louvain!

    ## 

    ## 6. Test OSLOM

    ## 

    ## Congratulation, you successfully install the version 2.5 of OSLOM!

The function has three parameters. `binpath` indicating the path to the
folder that will host the `bin` folder containing the binary files. By
default, the binary files are installed in R’s temporary directory
(`binpath = "tempdir"`). In this case the `bin` folder will be
automatically removed at the end of the R session. Alternatively, the
binary files can be installed in the bioregion’s package folder
(`binpath = "pkgfolder"`).

**Any other folder can be chosen but in any case PLEASE MAKE SURE to
update the `binpath` argument in
[netclu_infomap](https://bioRgeo.github.io/bioregion/reference/netclu_infomap.html),
[netclu_louvain](https://bioRgeo.github.io/bioregion/reference/netclu_louvain.html)
and
[netclu_oslom](https://bioRgeo.github.io/bioregion/reference/netclu_oslom.html)
accordingly.**

Note that in this case, and only in this case, `download_only` can be
set to `TRUE` (only the first step will be executed), but you must
ensure that the files have the necessary permissions to be executed as
programs otherwise.

The third parameter `infomap_version` indicating the Infomap version(s)
to install. Only the Infomap version 2.1.0, 2.6.0, 2.7.1 and 2.8.0 are
available for now.

The installation of the binary files is divided into six steps:

1.  Download bin.zip
2.  Unzip bin.zip in the binpath
3.  Check the permissions, try to change them automatically or propose
    an iterative process to change them manually
4.  Test that Infomap is running properly
5.  Test that Louvain is running properly
6.  Test that OSLOM is running properly

The function is designed to help you change the permissions during the
third step. Subsequently, the binary files of Infomap, Louvain, and
OSLOM are tested, and a file `check.txt` is added to each folder
`bin/INFOMAP/X.X.X`, `bin/LOUVAIN` and `bin/OSLOM` to inform the package
that
[netclu_infomap](https://bioRgeo.github.io/bioregion/reference/netclu_infomap.html),
[netclu_louvain](https://bioRgeo.github.io/bioregion/reference/netclu_louvain.html)
and
[netclu_oslom](https://bioRgeo.github.io/bioregion/reference/netclu_oslom.html)
can be used without any issues.

However, it may happen (particularly on macOS) that the function fails
to automatically or manually assist you in changing the permissions. In
this case, set the argument `download_only` to `TRUE` to download
`bin.zip` into a folder `myfolder` of your choice (different from
`tempdir` and `pkgfolder`) so that you can manually change the
permissions yourself. Note that you will be on your own; the package
will not be able to test the binary files before running
[netclu_infomap](https://bioRgeo.github.io/bioregion/reference/netclu_infomap.html),
[netclu_louvain](https://bioRgeo.github.io/bioregion/reference/netclu_louvain.html)
or
[netclu_oslom](https://bioRgeo.github.io/bioregion/reference/netclu_oslom.html)
(i.e., no `check.txt` will be added to the folders…). To bypass the
automatic check while using these functions, use the argument
`check_install = FALSE`, but without any guarantee of success.

``` r
install_binaries(binpath = "myfolder" ,
                 download_only =  TRUE,
                 infomap_version = c("2.1.0", 
                                     "2.6.0",
                                     "2.7.1",
                                     "2.8.0"))

comat <- matrix(sample(1000, 50), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

net <- similarity(comat, metric = "Simpson")
com <- netclu_infomap(net,
                      binpath = "myfolder",
                      check_install = TRUE)
```

## 2. Known issues

- The OpenMP versions of Infomap require `libomp-dev` on Ubuntu
  (`sudo apt-get install libomp-dev`) and `libomp` on macOS (install
  [Homebrew](https://brew.sh) and run `brew install libomp`).

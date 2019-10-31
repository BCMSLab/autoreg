# autoreg

A research compendium: Transcriptional regulation of autophagy during adipocyte differentiation

## Setting up the docker environment

The analysis was run on a [docker](https://hub.docker.com/r/bcmslab/autoreg/) image based on the the latest **rocker/verse**.
Other R packages were added to the image and were made available as an image that can be obtained and launched on any local
machine running [docker](https://hub.docker.com/r/bcmslab/autoreg/).

```bash
$ docker pull bcmslab/autoreg:latest
$ docker run -it bcmslab/autoreg:latest bash
```

## Obtaining the source code

The source code is hosted publicly on this repository in the form of a research compendium. This includes the scripts to
download, prepare and analyze the data. From within the container, [git](https://git-scm.com) can be
used to clone the source code.

The following code clones the repository containing the source code.

```bash
$ git clone http://github.com/BCMSLab/autoreg
```

## Runing the analysis

In the submodule `autoreg`, run `make`

```bash
$ cd autoreg
$ make analysis
```

## Details of the R environment
The version of **R** that was used to perform this analysis is the 3.7.0 (2019-04-07) on `x86\_64-pc-linux-gnu`.

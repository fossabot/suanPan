# <img src="Resource/suanPan-qr.svg" width="150" align="middle"/><img src="Resource/suanPan.svg" width="150" align="middle"/> suanPan

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1285221.svg)](https://doi.org/10.5281/zenodo.1285221)
[![License: GPL v3](https://img.shields.io/github/license/TLCFEM/suanPan.svg?color=44cc11)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/suanpan-manual/badge/?version=latest)](https://suanpan-manual.readthedocs.io/?badge=latest)
[![release](https://img.shields.io/github/release-pre/TLCFEM/suanPan.svg?color=44cc11)](https://github.com/TLCFEM/suanPan/releases)
[![suanpan](https://snapcraft.io//suanpan/badge.svg)](https://snapcraft.io/suanpan)
[![Chocolatey Version](https://img.shields.io/chocolatey/v/suanpan?color=44cc11)](https://chocolatey.org/packages/suanpan)
[![download](https://img.shields.io/github/downloads/TLCFEM/suanPan/total.svg?color=44cc11)](https://img.shields.io/github/downloads/TLCFEM/suanPan/total.svg?color=44cc11)
[![C/C++ CI](https://github.com/TLCFEM/suanPan/workflows/C/C++%20CI/badge.svg?branch=master)](https://github.com/TLCFEM/suanPan/actions)
[![AppVeyor branch](https://img.shields.io/appveyor/ci/TLCFEM/suanPan/master.svg?label=master&logo=appveyor)](https://ci.appveyor.com/project/TLCFEM/suanpan/branch/master)
[![Build Status](https://travis-ci.com/TLCFEM/suanPan.svg?branch=master)](https://travis-ci.com/TLCFEM/suanPan)
[![codecov](https://codecov.io/gh/TLCFEM/suanPan/branch/dev/graph/badge.svg)](https://codecov.io/gh/TLCFEM/suanPan)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/0754f66f6fae4829993eb879b222a32a)](https://www.codacy.com/app/TLCFEM/suanPan?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=TLCFEM/suanPan&amp;utm_campaign=Badge_Grade)
[![CodeFactor](https://www.codefactor.io/repository/github/tlcfem/suanpan/badge)](https://www.codefactor.io/repository/github/tlcfem/suanpan)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/TLCFEM/suanPan.svg?logo=lgtm)](https://lgtm.com/projects/g/TLCFEM/suanPan/context:cpp)
[![language](https://img.shields.io/github/languages/count/TLCFEM/suanPan.svg?color=44cc11)](https://github.com/TLCFEM/suanPan)
[![GitHub top language](https://img.shields.io/github/languages/top/TLCFEM/suanPan.svg?color=44cc11&logo=c%2B%2B)](https://github.com/TLCFEM/suanPan)
[![code-size](https://img.shields.io/github/languages/code-size/TLCFEM/suanPan.svg?color=44cc11)](https://img.shields.io/github/languages/code-size/TLCFEM/suanPan.svg?color=44cc11)
[![GitHub issues](https://img.shields.io/github/issues/TLCFEM/suanPan.svg?color=44cc11)](https://github.com/TLCFEM/suanPan/issues)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/TLCFEM/suanPan.svg?logo=lgtm)](https://lgtm.com/projects/g/TLCFEM/suanPan/alerts/)
[![Gitter](https://badges.gitter.im/suanPan-dev/community.svg)](https://gitter.im/suanPan-dev/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

## Introduction

[🧮 **suanPan**](https://tlcfem.github.io/suanPan/) is a finite element method (FEM) simulation platform for applications in fields such as solid mechanics and civil/structural/seismic engineering. The name **suanPan** (in some places such as suffix it is also abbreviated as **suPan**) comes from the term *Suan Pan* (算盤), which is [Chinese abacus](https://en.wikipedia.org/wiki/Suanpan). **suanPan** is written in modern high quality C++ code and is targeted to provide an efficient, concise, flexible and reliable FEM simulation platform.

**suanPan** is partially influenced by popular (non-)commercial FEA packages, such as [ABAQUS UNIFIED FEA](https://www.3ds.com/products-services/simulia/products/abaqus/), [ANSYS](http://www.ansys.com/) and [OpenSees](http://opensees.berkeley.edu/).

Please check documentation [here](https://tlcfem.gitbook.io/suanpan-manual/) and [here](http://suanpan-manual.rtfd.io/) for command references. Please consider star ⭐ the project!

## Features

The highlights of **suanPan** are

- **suanPan** is *fast*, both memory and thread safe.
- **suanPan** is designed based on the [shared memory](https://en.wikipedia.org/wiki/Shared_memory) model and supports parallelism on heterogeneous architectures, for example multi-threaded CPU + optional GPU. The parallelism is available for both element state updating and global matrix assembling.
- **suanPan** is open source and easy to be expanded to incorporate user-defined elements, materials, etc.
- **suanPan** separates the FEA model part from the linear algebra operation part, which significantly reduces the complexity of development.
- **suanPan** utilizes the new language features shipped with the latest standards (C++14, C++17, etc.), such as new STL containers, smart pointers and many others.
- **suanPan** supports simple visualization supported by [VTK](https://vtk.org/).

## Quick Start

Sample models are available for almost all models/commands. Please check the `Example` folder for details.

## Installation

All precompiled binaries can be executed out-of-box. Only 64-bit version is compiled. It is assumed that [AVX](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions) is available thus if the program fails, please check if your CPU supports AVX.

### Windows

The binaries, which are compiled with Intel MKL and VTK, are available on [Chocolatey](https://chocolatey.org/packages/suanpan), please use the following command to install the package.

1. Follow the [instructions](https://chocolatey.org/install) to install Chocolatey.

2. Make sure the [Microsoft Visual C++ Redistributable for Visual Studio 2019](https://aka.ms/vs/16/release/vc_redist.x64.exe) is installed. Alternatively, install it via Chocolatey.

    ```
    choco install vcredist140
    ```

3. Use the following command to install `suanPan`.

    ```
    choco install suanpan
    ```

### Linux

Linux users are strongly recommended to obtain the binaries via snap. The snap supports visualization via VTK and uses Intel MKL for linear algebra.

[![Get it from the Snap Store](https://snapcraft.io/static/images/badges/en/snap-store-black.svg)](https://snapcraft.io/suanpan)

[![asciicast](https://asciinema.org/a/341345.svg)](https://asciinema.org/a/341345)

### Other Platforms

Binaries precompiled with [OpenBLAS](https://www.openblas.net/) are provided by CI/CD on MacOS, Windows and Ubuntu. Please download the file from the [release](https://github.com/TLCFEM/suanPan/releases) page.

Advanced users can compile the program from source by themselves in order to enable

1. GPU based solvers which require available [MAGMA](http://icl.cs.utk.edu/magma/) library;
2. Visualization support provided by VTK library;
3. High performing linear algebra provided by Intel MKL library.

On Windows, to add file associations with `.sp` and `.supan` files, please run the `AddAssociation.bat` file with admin privilege. [Sublime Text](https://www.sublimetext.com/) autocompletion file is also provided. Please check the `Enhancement` folder.

## Dependency

Additional libraries used in **suanPan** are listed as follows.

- [**ARPACK**](https://www.caam.rice.edu/software/ARPACK/) version 0.96
- [**SPIKE**](http://www.ecs.umass.edu/~polizzi/spike/index.htm) version 1.0
- [**SuperLU**](https://portal.nersc.gov/project/sparse/superlu/) version 5.2.2 and [**SuperLU MT**](https://portal.nersc.gov/project/sparse/superlu/) version 3.1
- [**OpenBLAS**](https://github.com/xianyi/OpenBLAS) version 0.3.10
- [**TBB** Threading Building Blocks](https://github.com/oneapi-src/oneTBB) version 2020U2
- [**HDF5**](https://www.hdfgroup.org/solutions/hdf5/) version 1.10.6
- [**MUMPS**](http://mumps.enseeiht.fr/) version 5.2.1
- [**VTK**](https://vtk.org/) version 8.2
- [**MAGMA**](https://icl.cs.utk.edu/magma/) version 2.4.0
- [**Armadillo**](http://arma.sourceforge.net/) version 10.1
- [**Intel MKL**](https://software.intel.com/en-us/mkl)

Those libraries may depend on other libraries such as [zlib](https://zlib.net/) and [Szip](https://support.hdfgroup.org/doc_resource/SZIP/). Additional tools may be used by **suanPan**, they are

- [**UPX** the Ultimate Packer for eXecutables](https://upx.github.io/)

## How To Compile

Please refer to the corresponding [page](https://github.com/TLCFEM/suanPan-manual/blob/dev/docs/Tutorial/Compile.md) in manual for details.

## Happy Modelling

![an example of simulation of particle collision](Resource/particle-collision.gif)

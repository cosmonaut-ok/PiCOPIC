# PiCoPiC

*PiCoPiC* - PiCoPiC is a Code of Particle in Cell

[![Gitlab Pipeline Status](https://gitlab.com/my-funny-plasma/PIC/picopic/badges/master/pipeline.svg)](https://gitlab.com/my-funny-plasma/PIC/picopic/commits/master)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-green.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Language](https://img.shields.io/badge/Language-C%2B%2B%20%7C%20Python-green)
![Platform](https://img.shields.io/badge/platform-Linux-brightgreen)

## Motivation

*PiCoPiC* is a successor of *[PDP3](https://github.com/knuniv/pdp3)* (and its further [incarnation](https://github.com/cosmonaut-ok/pdp3)) project. Developed as electromagnetic particle-in-cell code for collisionless plasma simulation, PDP3 designed as single-thread code. Late versions of PDP3 are partially parallelized, but it can not be fully parallelized and can not be pupulated to large number of threads (it still effective with 2-8 threads). Such critical limitation of PDP3 was a motivation to replace the project as new fully parallel PIC code, called PiCoPiC.

## General Description
PiCoPiC still uses same algorithms as PDP3, but separates the simulation area to numerous domains, which are fully independent from each other, except synchronization points, when edges of their currents and fields are "stitched" with overlaying.

## Build

```bash
$ ./autogen.sh
$ ./configure [OPTIONS]
$ make
```

### Blackbox Tests

```bash
$ make test
```

## Configuration

*PiCoPiC* uses JSON format of configuration files (default configfile is `PiCoPiC.json`). Example of config located in the root of this repository.

## Output and data analysis

PiCoPiC follows the same principles for data output, as PDP3 data array output and python library and tools collection (python scripts and jupyter notebooks) for analysis of such data. Python library can be found in lib/python/picopic subdir and tools collection - in `tools` directory.

### System And Software Requirements

- C++ compiler: gcc 6.X or later (8.X recommended)
- autotools (autoconf, autoheader)
- git (to get sources)
- hdf5 library (1.10.5+ is recommmended)
- 'autoconf' and 'make' utils
- libgomp (OpenMP spec. version 3.0+).
- python3 + jupyter + matplotlib + numpy for visualization (anaconda - python scientific environment is recommended)
- doxygen, latex (texlive-full), imagemagick, pandoc, pandoc-citeproc _**(optional, for documentation only)**_

### HOWTO (linux, debian/ubuntu example)

#### 0. INSTALL PREREQUIRED SOFTWARE

**Install prerequired software with standard package management system**:

``` shell
root@host# apt-get install build-essential autoconf
```

**HDF5 library installation (optional)**:

If you going to use HDF5 format, you should build HDF5 library as prerequirement. As libhdf5 v.1.10.0 ships with debian, we need to download and compile fresh version (1.10.3) manually

<u>**WARNING!** HDF5 is default picopic output format. You must take care about HDF5 library presence (see below) or switch to plaintext format with `./configure` script (see "CONFIGURE AND BUILD PROJECT").</u>

* 0.a. install HDF5 with package management system (perform 0.a, 0.b or install in some other way, at your option):
``` shell
user@host$ sudo apt-get install libhdf5-cpp-100 libhdf5-dev
```

* 0.b. install HDF5 globally with "make" helper (perform project configuration first. See **CONFIGURE AND BUILD PROJECT** section):

``` shell
user@host$ sudo apt-get install wget
user@host$ sudo make hdf5-install
```

#### 1. **GIT CLONE PROJECT**

``` shell
user@host$ git clone https://gitlab.com/my-funny-plasma/PIC/picopic.git
user@host$ cd picopic
user@host$ git submodule update --init # require to enable external libraries
```

#### 2. **CONFIGURE AND BUILD PROJECT**

``` shell
# change your current directory to project's root directory
user@host$ cd /path/to/picopic/
user@host$ ./autogen.sh
user@host$ ./configure [OPTIONS] # use ./configure --help to view full options list
# 0.b.: optional step to install HDF5 library with built-in helper
user@host$ sudo make hdf5-install
# end of 0.b.
user@host$ make [COMPILE FLAGS] # see below about COMPILE FLAGS
```

<u> **WARNING!** HDF5 is default picopic output format. if you don't want to use HDF5 output format, you must disable it manually with `--disable-hdf5` option for `./configure` script </u>

#### 3. **TEST (optional)**

* Functional (end-to-end) testing:
```shell
user@host$ make test
```

* Unit testing (not implemented yet):
```shell
user@host$ make test-unit
```

#### 4. **RUN**

After compilation finished, you just need binary file `PiCoPiC` and `PiCoPiC.json`. You can copy this files to somewhere (if you skipped installation), edit `PiCoPiC.json` and run `./PiCoPiC`.

``` shell
user@host$ /path/to/PiCoPiC [ -h ] | [ -f /path/to/PiCoPiC.json ] # used PiCoPiC.json from current directory, if calling without any options
# or (much better), launch it as high-priority process
user@host$ nice -20 /path/to/PiCoPiC [ -h ] | [ -f /path/to/PiCoPiC.json ] # give some power to PiCoPiC!
```
> NOTE: it can take several days or weeks, and several hundred gigabytes of diskspace (yep, it is science, my deer friend).

#### 5. **VISUALIZATION (generate plots or animation)**

After your application finished modeling, you can build some visual model from generated data. Use python with matplotlib and numpy. [Anaconda](https://www.anaconda.com/download/#linux) as python distribution is recommended.

**...but... if you don't want to use anaconda...:**

``` shell
root@host# apt-get install python3 python3-matplotlib python3-numpy python3-scipy
```

**Pre-required software (for animation):**

``` shell
root@host# apt-get install ffmpeg
```

**Animation generation:**

``` shell
user@host$ /path/to/repository/with/picopic/tools/movie_3_E_RHObeam_r_z_map.py /path/to/data/directory [OPTIONS] # --help to get list of available options

```

**Interactive data analysis:**

Data analysis tools made as jupyter (ipython) notebooks. You can find them in `tools` subdir. Please, run jupyter from `tools` dir also

``` shell
user@host$ cd /path/to/PiCoPiC
user@host$ jupyter <notebook|lab> [tools/<notebook name>.ipynb]
```


**Non-interactive notebook launch:**

To run python notebook in non-interactive mode, please, use script `tools/nbrun.sh` from project root directory
``` shell
user@host$ /path/to/PiCoPiC/tools/nbrun.sh /path/to/PiCoPiC/tools/<notebook name`>.ipynb [data_path=\'/path/to/data/directory\'] [other_notebook_variable=other_notebook_value]
```

#### 6. **DOCUMENTATION GENERATION** (optional)

**Pre-required software:**

``` shell
root@host# doxygen texlive-latex-base texlive-latex-extra imagemagick graphviz
```
**Generate:**

``` shell
user@host$ make doc
```
Find documentation in:

- PDF
  `doc/app/latex/refman.pdf` - application
  `doc/vis/latex/refman.pdf` - visualization
- HTML
  `doc/app/html/index.xhtml` - application
  `doc/vis/html/index.xhtml` - visualization

### Hacking

#### Coding Style

- Indentication: 2 spaces
- Braces: BSD style
- Naming:
  - classes - CamelCase
  - files - camelCase, same as class name, begining from lower case
  - functions - snake_case, lower case
- Spacing:
  - after/around operators
  - no space after unary operators
  - no space before the postfix increment & decrement unary operators
  - no space around the . and -> structure member operators
- Declaring pointer data or a function that returns a pointer type: preferred to use of * is adjacent to the data name or function name and not adjacent to the type name
- Comments: '//' at line begining
- File naming: lowercase

#### GIT

##### Before working

1. Make sure, that you have clean repository and there are no unstaged untracked and uncommited files.

``` shell
user@host$ git status
```
> NOTE: If you have such files, you commit, remove or move out of your working directory (with PiCoPiC), than reset repository
> Also, you can just reset repository, or remove your local changes in other way
``` shell
user@host$ git reset --hard
```

2. Sync with upstream

``` shell
user@host$ git fetch origin
```

##### Working

1. Checkout to new branch

``` shell
user@host$ git checkout -b <your-branch-name> # use only `-` as word delimiters for better compatability
```

2. Make changes, write new code, add new features, experiment etc.

3. Commit your changes (do it as often, as you press `Ctrl+s` ;) ) and push to upstream

``` shell
user@host$ git add . # or git add your/selected/files
user@host$ git commit -m "<your comment>"
user@host$ git push origin <your-branch-name>
```

##### After working

When you finish some logical step of your work, you should merge your changes to master branch (as stable code). You can do it, using "Pull request" in github

1. Go to https://gitlab.com/my-funny-plasma/PIC/picopic/merge_requests

2. Press "New merge request" and choose "master" as base branch and "your-branch-name" as compare. Scroll to view changes against master branch

3. When all is ok, press "Create pull request" and confirm.

4. Wait for travis testing ends and (if they passed), merge your changes to master branch (squash and merge).

#### Debug

**Linux:**

1. build project with DEBUG option
``` shell
user@host$ ./autogen.sh && ./configure --enable-debug && make
```

2. run with gdb

``` shell
user@host$ gdb ./PiCoPiC
(gdb) run ## or perform some modifications first than run, e.g. set breakpoints
```
3. use experimental features

If you want to experiment and use unsafe features, disabled by default, use "--enable-experimental" key during configuration:

``` shell
user@host$ ./confugure --enable-experimental [other options] && make
```

#### Tools

Quick analytic calculator of plasma (aka Langmur) frequency, wake wavelength, Debye length etc. from PiCoPiC.json file
``` shell
user@host$ ./tools/quick_parameters_calculator.py <path/to/PiCoPiC.json>
```
See **VISUALIZATION** section also

You can edit jupyter notebooks with jupyter browser editor (opens with `jupyter notebook` or `jupyter lab` commands), atom (plugin: https://atom.io/packages/jupyter-notebook), emacs (plugin: https://github.com/millejoh/emacs-ipython-notebook), spyder (plugin: https://github.com/spyder-ide/spyder-notebook) etc

#### Bugs/Workarounds

> NOTE: It's recommended to disable HDF data file locking to avoid unreadable files, during incorrect application exit or other emergency stuff. Use `export HDF5_USE_FILE_LOCKING=FALSE` before launching PiCoPiC, jupyter or other tools to disable locking.

#### Known Issues

- some plasma instabilities with unknown nature found. Required to investigate root cause.

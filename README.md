# PicoPIC

*PicoPIC* - PicoPIC Is Code (Computer/Calculator) Of Particle In Cell

[![Gitlab Pipeline Status](https://gitlab.com/my-funny-plasma/PIC/picopic/badges/master/pipeline.svg)](https://gitlab.com/my-funny-plasma/PIC/picopic/commits/master)

## Motivation

*PicoPIC* is a successor of *[PDP3](https://github.com/cosmonaut-ok/pdp3/)* ([old repo](https://github.com/knuniv/pdp3)) project. Developed as electromagnetic particle-in-cell code for collisionless plasma simulation, PDP3 designed as single-thread code. Late versions of PDP3 are partially parallelized, but it can not be fully parallelized and can not be pupulated to large number of threads (it still effective with 2-8 threads). Such critical limitation of PDP3 was a motivation to replace the project as new fully parallel PIC code, called PicoPIC.

## General Description
PicoPIC still uses same algorithms as PDP3, but separates the simulation area to numerous areas, which are fully independent from each other, except synchronization points, when edges of their currents and fields are "stitched" with overlaying.

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

*PicoPIC* uses JSON format of configuration files (default configfile is `PicoPIC.json`). Example of config located in the root of this repository.

## Output and data analysis

PicoPIC follows the same principles for data output, as PDP3 data array output and python library and tools collection (python scripts and jupyter notebooks) for analysis of such data. Python library can be found in lib/python/picopic subdir and tools collection - in `tools` directory.


==================================================

# PDP3

### System And Software Requirements

- C++ compiler: gcc 6.X (or later, but not tested with 7.X and 8.X)
- git (to get sources)
- hdf5 library (1.10.1+ recommmended. 1.10.0 is buggy)
- 'autoconf' and 'make' utils
- libgomp (OpenMP spec. version 3.0+).
- python+jupyter+matplotlib+numpy for visualization (anaconda - python scientific environment is highly recommended)
- doxygen, latex (texlive-full), imagemagick, pandoc, pandoc-citeproc _**(optional, for documentation only)**_

### HOWTO (linux, debian/ubuntu example)

#### 0. INSTALL PREREQUIRED SOFTWARE

**Install prerequired software with standars package management system**:

``` shell
root@host# apt-get install build-essential autotools git
```

**HDF5 library installation (optional)**:

If you going to use HDF5 format, you should build HDF5 library as prerequirement. As libhdf5 v.1.10.0 ships with debian, we need to download and compile fresh version (1.10.3) manually

* 0.a. install HDF5 manually (perform 0.a, 0.b or install in some other way, at your option):
``` shell
user@host$ cd /tmp
user@host$ wget -qO- "https://www.hdfgroup.org/package/source-bzip/?wpdmdl=12594&refresh=5bbc7778635b21539078008" | tar xjf -
user@host$ cd /tmp/hdf5-1.10.3
user@host$ ./autogen.sh
user@host$ ./configure --enable-cxx --enable-build-mode=production --prefix=/usr/local
user@host$ sudo make install
user@host$ cd /tmp
user@host$ sudo rm -rf /tmp/hdf5-1.10.3
```

* 0.b. install HDF5 globally with make helper and include:\
You can install it with PDP3 build system's helper after application configuration, but before project compilation (see **COMPILE PROJECT** step)

#### 1. **GIT CLONE PROJECT**

``` shell
user@host$ git clone https://github.com/cosmonaut-ok/pdp3.git
user@host$ cd pdp3
user@host$ git submodule update --init # require to enable external libraries
```

#### 2. **CONFIGURE AND BUILD PROJECT**

``` shell
# change your current directory to project's root directory
user@host$ cd /path/to/pdp3/
user@host$ ./autogen.sh
user@host$ ./configure [--help|show help] [--some-options] # use ./configure --help to view full options list
# 0.b.: optional step to install HDF5 library with built-in helper
user@host$ make hdf5
user@host$ sudo make hdf5-install
# end of 0.b.
user@host$ make [COMPILE FLAGS] # see below about COMPILE FLAGS
```

#### 3. **TEST (optional)**

* Functional (end-to-end) testing:
```shell
user@host$ make test # or test-ext for extended testing (require more time)
```

* Unit testing:
```shell
user@host$ make test-unit
```

#### 4. **INSTALLATION (optional)**

You can install already compiled pdp3 with it's configfile (aka properties.xml) and required subdirs to separate directory. Just run:

``` shell
user@host$ make dist [RELEASE=/path/to/some/target/directory]
```

Than, you can go to this directory and run pdp3

#### 5. **RUN**

After compilation finished, you just need binary file `pdp3` and `parameters.xml`. You can copy this files to somewhere (if you skipped installation), edit `parameters.xml` and run `./pdp3`.

_**WARNING! pdp3 does not create data directory automatically. You need to create it before run (see `result_path` in parameters.xml)**_

``` shell
user@host$ mkdir pdp_result # or where you defined in parameters.xml. PDP3 does not use smth. like BOOST::filesystem to operate with directories
# NOTE: if you performed `make dist`, you just need `cd /path/to/target/dir`, edit `parameters.xml` (optional) and run `./pdp3`
user@host$ /path/to/pdp3 [ -h ] | [ -f /path/to/parameters.xml ] # used parameters.xml from current directory, if calling without any options
# or (much better), launch it as high-priority process
user@host$ nice -20 /path/to/pdp3 [ -h ] | [ -f /path/to/parameters.xml ] # give some power to pdp3!
```
> NOTE: it can take several days or weeks, and several hundred gigabytes of diskspace (yep, it is science, my deer friend).

#### 6. **VISUALIZATION (generate plots or animation)**

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
user@host$ /path/to/repository/with/pdp3/tools/movie_3_component.py /path/to/parameters.xml
# USE: /path/to/repository/with/pdp3/tools/movie_3_component.py -h to see list of all available options
```
> NOTE: you may use smth. like `python3 /path/to/repository/with/pdp3/tools/movie_3_component.py /path/to/parameters.xml`, if you don't use anaconda.

**Interactive data analysis:**

Data analysis tools made as jupyter (ipython) notebooks. You can find them in `tools` subdir. Please, run jupyter from `tools` dir also

``` shell
user@host$ cd /path/to/pdp3/tools
user@host$ jupyter <notebook|lab>
```

than you can choose jupyter notebook, corresponding to your needs and work with it.

**Non-interactive notebook launch:**

To run python notebook in non-interactive mode, please, use script `tools/nbrun.sh` from project root directory
``` shell
user@host$ /path/to/pdp3/tools/nbrun.sh /path/to/pdp3/tools/some.ipynb [config_path=\'/path/to/parameters.xml\'] [other_notebook_variable=other_notebook_value]
```

#### 7. **DOCUMENTATION GENERATION** (optional)

**Pre-required software:**

``` shell
root@host# doxygen texlive-latex-base texlive-latex-extra imagemagick
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
> NOTE: If you have such files, you commit, remove or move out of your working directory (with pdp3), than reset repository
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

1. Go to https://github.com/cosmonaut-ok/pdp3/pulls

2. Press "New pull request" and choose "master" as base branch and "your-branch-name" as compare. Scroll to view changes against master branch

3. When all is ok, press "Create pull request" and confirm.

4. Wait for travis testing ends and (if they passed), merge your changes to master branch (squash and merge).

#### Debug

**Linux:**

1. build project with DEBUG option
``` shell
user@host$ ./autogen.sh && ./configure --enable-debug && make
```

2. Set `debug` option to `true`,

3. Check out options
   - `particles->particles_kind->debug_number` (at least `1e4` is recommended)
   - `particles_bunch->debug_number` (at least `1e4` is recommended)
   - `geometry->debug_n_grid_r`
   - `geometry->debug_n_grid_z`
   - `file_save_parameters->debug_data_dump_interval`
in `parameters.xml` file before project run.

4. run with gdb

``` shell
user@host$ gdb ./pdp3
(gdb) run ## or perform some modifications first than run, e.g. set breakpoints
```
5. use experimental features

If you want to experiment and use unsafe features, disabled by default, use "-DEXPERIMENTAL" key in CFLAGS. For example:

``` shell
user@host$ ./confugure --enable-experimental && make
```

#### Tools

Quick analytic calculator of plasma (aka Langmur) frequency, wake wavelength, Debye length etc. from parameters.xml file
``` shell
user@host$ ./tools/quick_parameters_calculator.py <path/to/parameters.xml>
```
See **VISUALIZATION** section also

You can edit jupyter notebooks with jupyter browser editor (opens with `jupyter notebook` or `jupyter lab` commands), atom (plugin: https://atom.io/packages/jupyter-notebook), emacs (plugin: https://github.com/millejoh/emacs-ipython-notebook), spyder (plugin: https://github.com/spyder-ide/spyder-notebook) etc

#### Bugs/Workarounds

> NOTE: It's recommended to disable HDF data file locking to avoid unreadable files, during incorrect application exit or other emergency stuff. Use `export HDF5_USE_FILE_LOCKING=FALSE` before launching pdp3, jupyter or other tools to disable locking.

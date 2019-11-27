# PicoPIC

*PicoPIC* - PicoPIC Is Code (Computer/Calculator) Of Particle In Cell

## Motivation
*PicoPIC* is a successor of old *[PDP3](https://github.com/cosmonaut-ok/pdp3/)* project. Developed as electromagnetic particle-in-cell code for collisionless plasma simulation, PDP3 designed as single-thread code. Late versions of PDP3 are partially parallelized, but it can not be fully parallelized and can not be pupulated to large number of threads (it still effective with 2-8 threads). Such critical limitation of PDP3 was a motivation to replace the project as new fully parallel PIC code, called PicoPIC.

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
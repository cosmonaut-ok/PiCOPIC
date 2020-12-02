#!/bin/bash

rm -f data.h5 && make ${1} && ./PiCoPiC && ./tools/nbrun.sh ./tools/view_3_E_beam_density.ipynb data_path=\'/home/cosmonaut/dev/pic/picopic/data.h5\' timestamp=2e-12 clim_0=[-1e3,1e3] clim_1=[-1e3,1e3] clim_beam=[0,5e14]

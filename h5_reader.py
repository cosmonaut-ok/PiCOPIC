#!/usr/bin/python

import h5py
import numpy as np

f = h5py.File('/home/cosmonaut/dev/pic/picopic/data.h5', 'r')

dset = f['/E_r/rec/0-128_0-512/4']

shape = dset.shape
sum = np.reshape(dset[...], (shape[0] * shape[1]))

print(sum)

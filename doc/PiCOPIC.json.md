# Config file reference

### Generic

- bool: `"debug": true/false` - set or unset debugging output (not implemented yet. Only as compile option)
- bool: `"hdf5": true/false`, - enable or disable HDF5 output (not implemented yet. Only as compile option)
- float: `"macro_amount": 1e6` - set full amount of macroparticles (for all particle species and beams)
- integer: `"beam_plasma_macro_size_ratio": 2` - by default PiCOPIC makes that all of the macroparticles represents the same number of real particles. But, sometimes, lighter macroparticles for beams required (in this case larger number os such macroparticles should be present in the simulation area). This option allows to set ratio, between background plasma and beam macroparticles.

### geometry
- float `"radius": 0.075` - simulation area radius
- float `"longitude": 0.075` - simulation area longitude

#### grid_size
- integer `"radius": 128` - grid nodes amount by r-direction
- integer `"longitude": 128` - grid nodes amount by z-direction

#### areas_amount
Simulation area separated to several subareas (called "areas" also), which are independent and calculated in parallel. So, amount of such areas is the maximum number of parallel threads. Here can be defined amount of such areas by radius and longitude.

**INFO** areas amount automatically sets to 1 by r and 1 by z in case if application runs in single-thread mode.

- integer `"radius": 2`
- integer `"longitude": 4`

### pml

Set PML (doi:10.1006/jcph.1994.1159)

- float `"left_wall": 0`
- float `"right_wall": 0.05`
- float `"outer_wall": 0.2`
- float `"sigma_1": 1e-5`
- float `"sigma_2": 7e-2`

### time

- float `"start": 0` - simulation start time (0 as usual)
- float `"end": 3e-9` - simulation end time
- float `"step": 5e-14` - simulation time step

### particles

should be a list of particle species with options:
`"particles": [ { specie1 options ...}, { specie2 options ...}, ... ]`

**Options**:
- string `"name": "electrons"` - name of particles specie
- integer `"charge": -1` - particles charge (in electron charges)
- integer `"mass": 1` - particles mass (in electron masses)
- integer `"temperature": 1` - initial temperature (in electronvolts)

#### density

Set density of the particles specie ( $m^{-1} )

- float `"left": 1e+17,`
- float `"right": 1.05e+17`

### beams

should be a list of particle species with options:
`"beams": [ { beam1 options ...}, { beam1 options ...}, ... ]`

**Options**:
- string `"name": "electrons"` - name of the beam (used as beam_<name>)
- integer `"charge": -1` - charge of the beam (in electron charges)
- integer `"mass": 1` - mass of the beam (in electron masses)
- float `"velocity": 2.8e8` - initial velocity of the beam
- integer `"start_time": 0` - time, when beam injection begins

#### bunch

Single particles beam can contain multiple bunches with similar length and radius. Parameters of this bunches set here:

- integer `"amount": 2` - amount of bunches in the beam (use 1 for just one bunch) in meters
- float `"distance": 0.0914` - initial distnce between bunches (between end of leading and begin of following) in meters
- float `"length": 0.0056` - initial bunch length (beam length in case of single bunch) in meters
- float `"radius": 0.005` - bunch radius in meters
- float `"density": 1e+16` - initlal bunch density in meters


### data

- string `"data_root": "./simulation_result"` - directory, where simulation result dumps

#### compression

- some output engines (only HDF5 for now) supports on-the-fly compression. It useful for large data sets to decrease HDD load.

- bool `"use": true/false`
- integer `"level": 1` - compression level

#### probes

should be a list of particle species with options:
`"probes": [ { probe1 options ...}, { probe2 options ...}, ... ]`

- list of string-like options `"component": "temperature"` - component for probe to dump. Possible options are: **E\_r**, **E\_phi**, **E\_z**, **H\_r**, **H\_phi**, **H\_z**, **J\_r**, **J\_phi**, **J\_z**, **velocity** (only with specie), **position** (only with specie), **density** (only with specie), **temperature** (only with specie)
- string `"specie": "electrons"`, ("beam_<beam name>" for beams), optional, required only for density and temperature components
- integer "schedule": 20,
- list of string-like options `"shape": "dot"` - "rec" is for rectangle-shape probe, "row" - is for row-shaped probe "col" - is for column-shaped probe "dot" - is for dot-shaped probe
- integer `"r": 20` for dot, or "start->r" and "end->r" for rectangle and row
- integer `"z": 20` for dot, or "start->z" and "end->z" for rectangle and column

example:
```
 "start": { "r": 0, "z": 0 },
 "end": { "r": 128, "z": 512 }
```
### plot

Plot and Video sections used only for visualization part. It contains information about plot parameters. Currently it builds with python's matplotlib.

#### subplot

- string `"cmap": "terrain"` - colormap for all of the subplots (can be redefined in visualization scripts). Not implemented yet

##### ticks

Amount of the plot ticks by x, z axes and color bar

- integer `"x_amount": 20`
- integer `"y_amount": 8`
- integer `"z_amount": 20`
- integer `"cbar_amount": 3`

#### figure

- string `"color": "None"` - background color of the figure
- float `"width": 10.5` - width of the figure
- float `"height": 7` - height of the figure
- integer `"dpi": 100` - DPI of the figure

##### font

Figure and plots font parameters

- string `"family": "sans-serif"`
- integer `"size": 10`
- string `"name": "DejaVu Sans"`

#### video

Figure parameters

- string `"codec": "mjpeg"`
- integer `"fps": 20`
- integer `"bitrate": 32000`


## Typical example

```json
{
  "debug": false,
  "hdf5": false,
  "macro_amount": 2.1e5,
  "beam_plasma_macro_size_ratio": 2,
  "geometry": {
    "radius": 0.075,
    "longitude": 0.3,
    "grid_size": {
      "radius": 128,
      "longitude": 512
    },
    "areas_amount": {
      "radius": 2,
      "longitude": 8
    }
  },
  "pml": {
    "left_wall": 0,
    "right_wall": 0.05,
    "outer_wall": 0.2,
    "sigma_1": 1e-5,
    "sigma_2": 7e-2
  },
  "time": {
    "start": 0,
    "end": 3e-9,
    "step": 5e-14
  },
  "particles": [
    {
      "name": "electrons",
      "charge": -1,
      "mass": 1,
      "density": {
        "left": 1e+17,
        "right": 1.05e+17
      },
      "temperature": 1
    },
    {
      "name": "ions",
      "charge": 1,
      "mass": 1836,
      "density": {
        "left": 1e+17,
        "right": 1.05e+17
      },
      "temperature": 0.1
    }
  ],
  "beams": [
    {
      "name": "electrons",
      "charge": -1,
      "mass": 1,
      "velocity": 2.8e8,
      "start_time": 0,
      "bunch": {
        "amount": 2,
        "distance": 0.0914,
        "length": 0.0056,
        "radius": 0.005,
        "density": 1e+16
      }
    }
  ],
  "data": {
    "data_root": "./simulation_result",
    "compression": {
      "use": false,
      "level": 1
    },
    "probes": [
      {
        "component": "temperature",
        "specie": "electrons",
        "schedule": 20,
          "shape": "dot",
	  "r": 20,
	  "z": 40
      },
      {
        "component": "E_r",
        "schedule": 20,
        "shape": "rec",
        "start": { "r": 0, "z": 0 },
        "end": { "r": 128, "z": 512 }
      },
      {
        "component": "E_phi",
        "schedule": 20,
        "shape": "rec",
        "start": { "r": 0, "z": 0 },
        "end": { "r": 128, "z": 512 }
      },
      {
        "component": "E_z",
        "schedule": 20,
        "shape": "rec",
        "start": { "r": 0, "z": 0 },
        "end": { "r": 128, "z": 512 }
      },
      {
        "component": "H_r",
        "schedule": 20,
        "shape": "rec",
        "start": { "r": 0, "z": 0 },
        "end": { "r": 128, "z": 512 }
      },
      {
        "component": "H_phi",
        "schedule": 20,
        "shape": "rec",
        "start": { "r": 0, "z": 0 },
        "end": { "r": 128, "z": 512 }
      },
      {
        "component": "H_z",
        "schedule": 20,
        "shape": "rec",
        "start": { "r": 0, "z": 0 },
        "end": { "r": 128, "z": 512 }
      },
      {
        "component": "J_r",
        "schedule": 20,
        "shape": "rec",
        "start": { "r": 0, "z": 0 },
        "end": { "r": 128, "z": 512 }
      },
      {
        "component": "J_phi",
        "schedule": 20,
        "shape": "rec",
        "start": { "r": 0, "z": 0 },
        "end": { "r": 128, "z": 512 }
      },
      {
        "component": "J_z",
        "schedule": 20,
        "shape": "rec",
        "start": { "r": 0, "z": 0 },
        "end": { "r": 128, "z": 512 }
      },
      {
        "component": "velocity",
        "specie": "electrons",
        "schedule": 10,
        "shape": "rec",
        "start": { "r": 0, "z": 0 },
        "end": { "r": 128, "z": 512 }
      },
      {
        "component": "position",
        "specie": "electrons",
        "schedule": 10,
        "shape": "rec",
        "start": { "r": 0, "z": 0 },
        "end": { "r": 128, "z": 512 }
      },
      {
        "component": "velocity",
        "specie": "beam_electrons",
        "schedule": 10,
        "shape": "rec",
        "start": { "r": 0, "z": 0 },
        "end": { "r": 128, "z": 512 }
      },
      {
        "component": "position",
        "specie": "beam_electrons",
        "schedule": 10,
        "shape": "rec",
        "start": { "r": 0, "z": 0 },
        "end": { "r": 128, "z": 512 }
      }
    ],
    "plot": {
      "subplot": {
        "cmap": "terrain",
        "ticks": {
          "x_amount": 20,
          "y_amount": 8,
          "z_amount": 20,
          "cbar_amount": 3
        }
      },
      "figure": {
        "color": "None",
        "width": 10.5,
        "height": 7,
        "dpi": 100,
        "font": {
          "family": "sans-serif",
          "size": 10,
          "name": "DejaVu Sans"
        }
      },
      "video": {
        "codec": "mjpeg",
        "fps": 20,
        "bitrate": 32000
      }
    }
  }
}
```

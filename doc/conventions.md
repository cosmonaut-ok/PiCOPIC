# Internal conventions

## Naming

1. Files: camelcase with small 1st letter
2. classes: camelcase with capital 1st letter
3. functions/methods: snake case without capital letters

## Defines naming

1. `ENABLE_`-convension: set define with `ENABLE_` prefix for on/off togglers
2. `_OPTION`-convension: set internal variable in `configure.ac` script bound to  `ENABLE_` define with the same name and `_OPTION` suffix.
3. `SWITCH_`-convension: set define with `SWITCH_` prefix to define from some variants (e.g. `#define SWITCH_PUSHER_BORIS`
4. Lockers: set of development and experimental options should be protected with `--enable-experimental` and `--enable-development` configure flags.


## Other

1. Calculation with writing to class member variables should be performed inside of this class (e.g. positions and velocities - Particles class, currents - Current class). Exception is AGlue, which performs cross-domains interaction.
2. Parts of complex names should be follow from general to specific (e.g. FieldElectric for electric field).
3. Never set `using namespace xxxx` in header files. It can cause unpredicted behavior and namespace conflicts at inclusion


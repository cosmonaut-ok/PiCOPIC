# Internal conventions

## Naming

1. Files: camelcase with small 1st letter
2. classes: camelcase with capital 1st letter
3. functions/methods: snake case without capital letters

## Other

1. Calculation with writing to class member variables should be performed inside of this class (e.g. positions and velocities - Particles class, currents - Current class). Exception is AGlue, which performs cross-domains interaction.
2. Parts of complex names should be follow from general to specific (e.g. FieldElectric for electric field).
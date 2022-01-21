# External fields configuration

PiCoPiC.json should contain `external_fields` section. with `electric` and `magnetic` subsections. Each of them should contain `formfactor` and `params` option. `formfactor` should be one of selections (e.g. `const`, `const_r`, `const_z`, `linegrad`, `linegrad_r`, `linegrad_z` etc) and `params` should contain array of params, specific for each formfactor (see. the table of formfactors and its parameters below).

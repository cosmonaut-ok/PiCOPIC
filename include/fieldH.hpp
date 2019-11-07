#pragma once

#include "field.hpp"
#include "fieldE.hpp"
#include "specieP.hpp"

class SpecieP; // define dummy class to avoid an error while loading real links
class FieldE;

class FieldH : public Field
{
public:
  Grid3D<double> field_at_et;
  vector<SpecieP *> species_p;
  FieldE* field_e;

public:
  FieldH(Geometry *geom, TimeSim *t, vector<SpecieP *> species);
  FieldH() {};
  ~FieldH(void) {};

  void calc_field_cylindrical();
  vector3d<double> get_field(double radius, double longitude);

  // HField(Geometry *geom1);
  // HField(void);
  // ~HField(void);
  // void calc_field(EField *e_field1, Time *time1);
  // void set_homogeneous_h(double E_r, double E_phi, double E_z);
  // double* get_field(double radius, double longitude);

  // double *get_1d_field_r();
  // double *get_1d_field_phi();
  // double *get_1d_field_z();
};

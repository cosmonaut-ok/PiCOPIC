#ifndef TINYVEC3D_H_INCLUDED
#define TINYVEC3D_H_INCLUDED
#endif

#include <sstream>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <cmath>

extern "C"
{
  namespace tinyvec3d
  {
    double* mkvector3d(double x0, double x1, double x2);
    double* tv_copy(double* v_this);
    void tv_copy_components(double* v_this, double* v_other);
    void tv_cross(double* v_this, double* v_other);
    double tv_dot(double* v_this, double* v_other);
    void tv_add(double* v_this, double* v_other);
    void tv_subst(double* v_this, double* v_other);
    void tv_plus(double* v_this, double number);
    void tv_minus(double* v_this, double number);
    void tv_product(double* v_this, double number);
    void tv_div(double* v_this, double number);
    void tv_power(double* v_this, double number);
    void tv_abs(double* v_this);
    void tv_abs2(double* v_this);
    double tv_squared_sum(double* v_this);
  }
}

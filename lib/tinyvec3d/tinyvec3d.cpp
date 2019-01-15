#include "tinyvec3d.h"
//! tiny library, which implements operations
//! with 3-dimention mathematical vectors
//! Library built as wrapper over c-style
//! 3-component arrays [ x1, x2, x3 ]
//! and does not contain any protections
//! from incorrect length or other incorrect
//! data. Because it is TINY :) Enjoy.
//! Be patient almost operations are destructive
//! and affects first given vector

#include "tinyvec3d.h"

namespace tinyvec3d {

  double* mkvector3d(double x0, double x1, double x2)
  {
    double* vector = new double[3];
    vector[0] = x0;
    vector[1] = x1;
    vector[2] = x2;

    return vector;
  }

  double* tv_copy(double* v_this)
  {
    double* newvec = new double[3];
    for (int i=0; i < 3; i++)
      newvec[i] = v_this[i];

    return newvec;
  }

  void tv_copy_components(double* v_this, double* v_other)
  {
    for (int i=0; i < 3; i++)
      v_other[i] = v_this[i];
  }

  void tv_cross(double* v_this, double* v_other)
  {
    double* vtmp = tv_copy(v_this);

    vtmp[0] = v_this[1] * v_other[2] - v_this[2] * v_other[1];
    vtmp[1] = v_this[2] * v_other[0] - v_this[0] * v_other[2];
    vtmp[2] = v_this[0] * v_other[1] - v_this[1] * v_other[0];

    for (int i=0; i < 3; i++)
      v_this[i] = vtmp[i];

    delete [] vtmp;
  }

  double tv_dot(double* v_this, double* v_other)
  {
    double res = 0;
    for (int i=0; i < 3; i++)
      res += v_this[i] * v_other[i];

    return res;
  }

  void tv_add(double* v_this, double* v_other)
  {
    for (int i=0; i < 3; i++)
      v_this[i] += v_other[i];
  }

  void tv_subst(double* v_this, double* v_other)
  {
    for (int i=0; i < 3; i++)
      v_this[i] -= v_other[i];
  }
  // TODO: vector-number: product, div, power, abs, abs2, size_t, floor

  ////////////////////////////////////////////////////////////////////////////////

  // vector-number
  void tv_plus(double* v_this, double number)
  {
    for (int i=0; i < 3; i++)
      v_this[i] += number;
  }

  void tv_minus(double* v_this, double number)
  {
    for (int i=0; i < 3; i++)
      v_this[i] -= number;
  }

  void tv_product(double* v_this, double number)
  {
    for (int i=0; i < 3; i++)
      v_this[i] *= number;
  }

  void tv_div(double* v_this, double number)
  {
    for (int i=0; i < 3; i++)
      v_this[i] /= number;
  }

  void tv_power(double* v_this, double number)
  {
    for (int i=0; i < 3; i++)
      v_this[i] = std::pow(v_this[i], number);
  }

  void tv_abs(double* v_this)
  {
    for (int i=0; i < 3; i++)
      v_this[i] = std::abs(v_this[i]);
  }

  void tv_abs2(double* v_this)
  {
    tv_power(v_this, 2);
  }

  double tv_squared_sum(double* v_this)
  {
    return tv_dot(v_this, v_this);
  }

  // TODO: vector-vector: max, min,
}

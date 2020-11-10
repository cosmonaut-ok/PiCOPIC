/* 
 * This file is part of the PiCoPiC distribution (https://github.com/cosmonaut-ok/PiCoPiC).
 * Copyright (c) 2020 Alexander Vynnyk.
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _VECTOR_3D_HPP_
#define _VECTOR_3D_HPP_

#include <math.h>
#include <vector>

#include "defines.hpp"
#include "msg.hpp"

template <typename T>
class vector3d
{
  std::vector<T> vec3d{0, 0, 0};

public:
  vector3d() // create zero-vector, if other is not given
  {
    vec3d[0] = 0;
    vec3d[1] = 0;
    vec3d[2] = 0;
  };

  vector3d(T x0, T x1, T x2)
  {
    vec3d[0] = x0;
    vec3d[1] = x1;
    vec3d[2] = x2;
  };

  // add operators
  // call vector3d components as array elements
  T& operator[] (int index)
  {
    return vec3d[index];
  };

  // == operator
  bool operator== (vector3d const& rhs)
  {
    return (vec3d[0] == rhs.vec3d[0]
            && vec3d[1] == rhs.vec3d[1]
            && vec3d[2] == rhs.vec3d[2]);
  };

  bool operator== (T const& rhs)
  {
    return (vec3d[0] == rhs
            && vec3d[1] == rhs
            && vec3d[2] == rhs);
  };

  // != operator
  bool operator!= (vector3d const& rhs)
  {
    return (vec3d[0] != rhs.vec3d[0]
            || vec3d[1] != rhs.vec3d[1]
            || vec3d[2] != rhs.vec3d[2]);
  };

  bool operator!= (T const& rhs)
  {
    return (vec3d[0] != rhs
            || vec3d[1] != rhs
            || vec3d[2] != rhs);
  };

  // = for number and another vector3d
  vector3d& operator= (T const& rhs)&
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      vec3d[i] = rhs;

    return *this;
  };

  vector3d& operator= (vector3d const& rhs)&
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      vec3d[i] = rhs.vec3d[i];

    return *this;
  };

  // += for number and another vector3d
  vector3d& operator+= (T const& rhs)&
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      vec3d[i] += rhs;

    return *this;
  };

  vector3d& operator+= (vector3d const& rhs)&
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      vec3d[i] += rhs.vec3d[i];

    return *this;
  };

  // -= for number and another vector3d
  vector3d& operator-= (T const& rhs)&
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      vec3d[i] -= rhs;

    return *this;
  };

  vector3d& operator-= (vector3d const& rhs)&
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      vec3d[i] -= rhs.vec3d[i];

    return *this;
  };

  // *= for number and another vector3d
  vector3d& operator*= (T const& rhs)&
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      vec3d[i] *= rhs;

    return *this;
  };

  vector3d& operator*= (vector3d const& rhs)&
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      vec3d[i] *= rhs.vec3d[i];

    return *this;
  };

  // /= for number and another vector3d
  vector3d& operator/= (T const& rhs)&
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      vec3d[i] /= rhs;

    return *this;
  };

  vector3d& operator/= (vector3d const& rhs)&
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      vec3d[i] /= rhs.vec3d[i];

    return *this;
  };

  friend vector3d operator+ (vector3d lhs, const vector3d& rhs)
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      lhs.vec3d[i] += rhs.vec3d[i];

    return lhs;
  };

  friend vector3d operator+ (vector3d lhs, const T& rhs)
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      lhs.vec3d[i] += rhs;

    return lhs;
  };

  friend vector3d operator- (vector3d lhs, const vector3d& rhs)
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      lhs.vec3d[i] -= rhs.vec3d[i];

    return lhs;
  };

  friend vector3d operator- (vector3d lhs, const T& rhs)
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      lhs.vec3d[i] -= rhs;

    return lhs;
  };

  friend vector3d operator* (vector3d lhs, const vector3d& rhs)
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      lhs.vec3d[i] *= rhs.vec3d[i];

    return lhs;
  };

  friend vector3d operator* (vector3d lhs, const T& rhs)
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      lhs.vec3d[i] *= rhs;

    return lhs;
  };

  friend vector3d operator/ (vector3d lhs, const vector3d& rhs)
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      lhs.vec3d[i] /= rhs.vec3d[i];

    return lhs;
  };

  friend vector3d operator/ (vector3d lhs, const T& rhs)
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      lhs.vec3d[i] /= rhs;

    return lhs;
  };

  friend vector3d operator% (vector3d lhs, const vector3d& rhs)
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      lhs.vec3d[i] %= rhs.vec3d[i];

    return lhs;
  };

  friend vector3d operator% (vector3d lhs, const T& rhs)
  {
#pragma omp simd
    for (unsigned short i=0; i<3; ++i)
      lhs.vec3d[i] %= rhs;

    return lhs;
  };

  // methods
  T dot(const vector3d rhs)
  {
    T ret = 0;
#pragma omp simd
    for (int i = 0; i < 3; i++)
      ret += vec3d[i] * rhs.vec3d[i];

    return T(ret);
  };

  vector3d cross(vector3d rhs)
  {
    vector3d vtmp = *this;

    vtmp[0] = vec3d[1] * rhs.vec3d[2] - vec3d[2] * rhs.vec3d[1];
    vtmp[1] = vec3d[2] * rhs.vec3d[0] - vec3d[0] * rhs.vec3d[2];
    vtmp[2] = vec3d[0] * rhs.vec3d[1] - vec3d[1] * rhs.vec3d[0];

    for (int i=0; i < 3; i++)
      vec3d[i] = vtmp[i];

    return *this;
  };

  T length2()
  {
    return T(vec3d[0]*vec3d[0] + vec3d[1]*vec3d[1] + vec3d[2]*vec3d[2] );
  };

  T length()
  {
    return T(sqrt(length2()));
  };

  vector3d pow(T power)
  {
#pragma omp simd
    for (int i=0; i < 3; i++)
      vec3d[i] = std::pow(vec3d[i], power);

    return *this;
  };

  vector3d abs()
  {
#pragma omp simd
    for (int i=0; i < 3; i++)
      vec3d[i] = std::abs(vec3d[i]);

    return *this;
  };
};

#endif // end of _VECTOR_3D_HPP_

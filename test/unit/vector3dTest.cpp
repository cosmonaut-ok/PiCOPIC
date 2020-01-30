#include <gtest/gtest.h>
#include "math/vector3d.hpp"

namespace {
  TEST(vector3d, constructor)
  {
    vector3d<int> vec(1, -2, 3);

    EXPECT_EQ(vec[0], 1);
    EXPECT_EQ(vec[1], -2);
    EXPECT_EQ(vec[2], 3);
  }

  TEST(vector3d, assign)
  {
    vector3d<int> vec(1, -2, 3);
    vector3d<int> rhs(-1, 2, -3);
    int a = 10;
    double b = 10.1;

    vec = a;
    EXPECT_EQ(vec[0], 10);
    EXPECT_EQ(vec[1], 10);
    EXPECT_EQ(vec[2], 10);

    vec = b;
    EXPECT_NE(vec[0], 10.1);
    EXPECT_NE(vec[1], 10.1);
    EXPECT_NE(vec[2], 10.1);

    vec = rhs;
    EXPECT_EQ(vec[0], -1);
    EXPECT_EQ(vec[1], 2);
    EXPECT_EQ(vec[2], -3);
  }

  TEST(vector3d, equality)
  {
    vector3d<int> vec(1, -2, 3);
    vector3d<int> rhs(1, -2, 3);
    vector3d<int> rhs_fail(-1, -2, -3);

    // positive
    EXPECT_TRUE(vec == rhs);
    EXPECT_FALSE(vec == rhs_fail);

    // not equal
    EXPECT_FALSE(vec != rhs);
    EXPECT_TRUE(vec != rhs_fail);

    // positive with number
    int a = 10;
    vec = a;
    EXPECT_FALSE(vec[0] != 10);

    // negative with number
    int b = 11;
    vec = b;
    EXPECT_TRUE(vec[0] != 10);
  }


  TEST(vector3d, add_assign)
  {
    vector3d<int> vec0(1, -2, 3);
    vector3d<int> vec1(1, -2, 3);
    vector3d<int> vec2(1, -2, 3);
    vector3d<int> vec3(1, -2, 3);
    vector3d<int> rhs(1, 2, -2);
    vector3d<int> rhs_fail(-1, -2, -3);
    int a = 10;
    int b = -10;
    double c = 10.1;

    vec0 += rhs;
    EXPECT_EQ(vec0[0], 2);
    EXPECT_EQ(vec0[1], 0);
    EXPECT_EQ(vec0[2], 1);

    vec1 += a;
    EXPECT_EQ(vec1[0], 11);
    EXPECT_EQ(vec1[1], 8);
    EXPECT_EQ(vec1[2], 13);

    vec2 += b;
    EXPECT_EQ(vec2[0], -9);
    EXPECT_EQ(vec2[1], -12);
    EXPECT_EQ(vec2[2], -7);

    vec3 += c;
    EXPECT_EQ(vec3[0], 11);
    EXPECT_EQ(vec3[1], 8);
    EXPECT_EQ(vec3[2], 13);
  }

  TEST(vector3d, substract_assign)
  {
    vector3d<int> vec0(1, -2, 3);
    vector3d<int> vec1(1, -2, 3);
    vector3d<int> vec2(1, -2, 3);
    vector3d<int> vec3(1, -2, 3);
    vector3d<int> rhs(1, 2, -2);
    vector3d<int> rhs_fail(-1, -2, -3);
    int a = 10;
    int b = -10;
    double c = 10.1;

    vec0 -= rhs;
    EXPECT_EQ(vec0[0], 0);
    EXPECT_EQ(vec0[1], -4);
    EXPECT_EQ(vec0[2], 5);

    vec1 -= a;
    EXPECT_EQ(vec1[0], -9);
    EXPECT_EQ(vec1[1], -12);
    EXPECT_EQ(vec1[2], -7);

    vec2 -= b;
    EXPECT_EQ(vec2[0], 11);
    EXPECT_EQ(vec2[1], 8);
    EXPECT_EQ(vec2[2], 13);

    vec3 -= c;
    EXPECT_EQ(vec3[0], -9);
    EXPECT_EQ(vec3[1], -12);
    EXPECT_EQ(vec3[2], -7);
  }

  TEST(vector3d, multiply_assign)
  {
    vector3d<int> vec0(1, -2, 3);
    vector3d<int> vec1(1, -2, 3);
    vector3d<int> vec2(1, -2, 3);
    vector3d<int> vec3(1, -2, 3);
    vector3d<int> rhs(1, 2, -2);
    vector3d<int> rhs_fail(-1, -2, -3);
    int a = 10;
    int b = -10;
    double c = 10.1;

    vec0 *= rhs;
    EXPECT_EQ(vec0[0], 1);
    EXPECT_EQ(vec0[1], -4);
    EXPECT_EQ(vec0[2], -6);

    vec1 *= a;
    EXPECT_EQ(vec1[0], 10);
    EXPECT_EQ(vec1[1], -20);
    EXPECT_EQ(vec1[2], 30);

    vec2 *= b;
    EXPECT_EQ(vec2[0], -10);
    EXPECT_EQ(vec2[1], 20);
    EXPECT_EQ(vec2[2], -30);

    vec3 *= c;
    EXPECT_EQ(vec3[0], 10);
    EXPECT_EQ(vec3[1], -20);
    EXPECT_EQ(vec3[2], 30);
  }

  TEST(vector3d, division_assign)
  {
    vector3d<int> vec0(1, -2, 3);
    vector3d<int> vec1(1, -2, 3);
    vector3d<int> vec2(1, -2, 3);
    vector3d<int> vec3(1, -2, 3);
    vector3d<int> rhs(1, 2, -2);
    vector3d<int> rhs_fail(-1, -2, -3);
    int a = 10;
    int b = -10;
    double c = 10.1;

    vec0 /= rhs;
    EXPECT_EQ(vec0[0], 1);
    EXPECT_EQ(vec0[1], -1);
    EXPECT_EQ(vec0[2], -1);

    vec1 /= a;
    EXPECT_EQ(vec1[0], 0);
    EXPECT_EQ(vec1[1], 0);
    EXPECT_EQ(vec1[2], 0);

    vec2 /= b;
    EXPECT_EQ(vec2[0], 0);
    EXPECT_EQ(vec2[1], 0);
    EXPECT_EQ(vec2[2], 0);

    vec3 /= c;
    EXPECT_EQ(vec3[0], 0);
    EXPECT_EQ(vec3[1], 0);
    EXPECT_EQ(vec3[2], 0);
  }

  TEST(vector3d, dot)
  {
    vector3d<int> vec(1, -2, 3);
    vector3d<int> vec2(-2, 3, 1);

    int dot = vec.dot(vec2);

    EXPECT_EQ(dot, -5);
  }


  TEST(vector3d, cross)
  {
    vector3d<int> vec(1, -2, 3);
    vector3d<int> vec2(-2, 3, 1);

    vector3d<int> cross = vec.cross(vec2);

    EXPECT_EQ(cross[0], -11);
    EXPECT_EQ(cross[1], -7);
    EXPECT_EQ(cross[2], -1);
  }

  TEST(vector3d, squared_length)
  {
    vector3d<int> vec(1, -2, 3);

    int length2 = vec.length2();

    EXPECT_EQ(length2, 14);
  }

  TEST(vector3d, power)
  {
    vector3d<int> vec(1, -2, 3);
    int power = 10;

    vector3d<int> pow = vec.pow(power);

    EXPECT_EQ(pow[0], 1);
    EXPECT_EQ(pow[1], 1024);
    EXPECT_EQ(pow[2], 59049);
  }

  TEST(vector3d, abs)
  {
    vector3d<int> vec(1, -2, 3);

    vector3d<int> abs = vec.abs();

    EXPECT_EQ(abs[0], 1);
    EXPECT_EQ(abs[1], 2);
    EXPECT_EQ(abs[2], 3);
  }
}

// SPDX-License-Identifier: MIT

/// \file
/// \brief Unit tests for dfe::poly

#include <array>
#include <valarray>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "dfe/dfe_poly.hpp"

// make std::pair printable to allow direct comparision w/ boost::test

namespace std {
template<typename T>
inline ostream&
operator<<(ostream& os, const pair<T, T>& x) {
  os << "{" << x.first << ", " << x.second << "}";
  return os;
}
} // namespace std

// test different input container types

#define COEFFS \
  { 1.0, 2.0, 0.25, 0.025 }
constexpr double X0 = 0.5;
constexpr double Y0 = 2.065625;
constexpr double D0 = 2.2687500000000003;
constexpr std::pair<double, double> YD0 = {Y0, D0};

BOOST_AUTO_TEST_CASE(poly_initializerlist) {
  BOOST_TEST(dfe::polynomial_val(X0, COEFFS) == Y0);
  BOOST_TEST(dfe::polynomial_der(X0, COEFFS) == D0);
  BOOST_TEST(dfe::polynomial_valder(X0, COEFFS) == YD0);
}

BOOST_AUTO_TEST_CASE(poly_array) {
  double coeffs[] = COEFFS;
  BOOST_TEST(dfe::polynomial_val(X0, coeffs) == Y0);
  BOOST_TEST(dfe::polynomial_der(X0, coeffs) == D0);
  BOOST_TEST(dfe::polynomial_valder(X0, coeffs) == YD0);
}

BOOST_AUTO_TEST_CASE(poly_stdarray) {
  std::array<double, 4> coeffs = COEFFS;
  BOOST_TEST(dfe::polynomial_val(X0, coeffs) == Y0);
  BOOST_TEST(dfe::polynomial_der(X0, coeffs) == D0);
  BOOST_TEST(dfe::polynomial_valder(X0, coeffs) == YD0);
}

BOOST_AUTO_TEST_CASE(poly_stdvector) {
  std::vector<double> coeffs = COEFFS;
  BOOST_TEST(dfe::polynomial_val(X0, coeffs) == Y0);
  BOOST_TEST(dfe::polynomial_der(X0, coeffs) == D0);
  BOOST_TEST(dfe::polynomial_valder(X0, coeffs) == YD0);
}

// use std::valarray to calculate polynomial for multiple x values at once

BOOST_AUTO_TEST_CASE(poly_valarray) {
  std::valarray<float> x(1024);
  for (std::size_t i = 0; i < x.size(); ++i) {
    x[i] = -1.0 + (2.0 / x.size()) * i;
  }
  auto linear = dfe::polynomial_val(x, {0.0, 1.0});
  auto quadratic = dfe::polynomial_val(x, {0.5, 0.0, 1.0});

  BOOST_TEST(x.size() == linear.size());
  BOOST_TEST(x.size() == quadratic.size());
  BOOST_TEST((x - linear).sum() == 0.0);
  BOOST_TEST((0.5 + x * x - quadratic).sum() == 0.0);
}

// special case of empty coefficients

BOOST_AUTO_TEST_CASE(poly_empty) {
  BOOST_TEST(dfe::polynomial_val(-1.0, std::array<double, 0>{}) == 0.0);
  BOOST_TEST(dfe::polynomial_val(+0.0, std::array<double, 0>{}) == 0.0);
  BOOST_TEST(dfe::polynomial_val(+1.0, std::array<double, 0>{}) == 0.0);
}

// test different fixed polynomial orders

BOOST_AUTO_TEST_CASE(poly_const) {
  BOOST_TEST(dfe::polynomial_val(-1.0, {42.0}) == 42.0);
  BOOST_TEST(dfe::polynomial_val(+0.0, {42.0}) == 42.0);
  BOOST_TEST(dfe::polynomial_val(+1.0, {42.0}) == 42.0);
}

BOOST_AUTO_TEST_CASE(poly_linear) {
  BOOST_TEST(dfe::polynomial_val(-0.5, {42.0, 1.0}) == 41.5);
  BOOST_TEST(dfe::polynomial_val(+0.0, {42.0, 1.0}) == 42.0);
  BOOST_TEST(dfe::polynomial_val(+0.5, {42.0, 1.0}) == 42.5);
}

BOOST_AUTO_TEST_CASE(poly_quadratic) {
  BOOST_TEST(dfe::polynomial_val(-0.5, {42.0, 1.0, 0.5}) == 41.625);
  BOOST_TEST(dfe::polynomial_val(+0.0, {42.0, 1.0, 0.5}) == 42.0);
  BOOST_TEST(dfe::polynomial_val(+0.5, {42.0, 1.0, 0.5}) == 42.625);
}

BOOST_AUTO_TEST_CASE(poly_cubic) {
  BOOST_TEST(dfe::polynomial_val(-0.5, {42.0, 1.0, 0.5, -1.0}) == 41.75);
  BOOST_TEST(dfe::polynomial_val(+0.0, {42.0, 1.0, 0.5, -1.0}) == 42.0);
  BOOST_TEST(dfe::polynomial_val(+0.5, {42.0, 1.0, 0.5, -1.0}) == 42.5);
}

// SPDX-License-Identifier: MIT
// Copyright 2018-2019 Moritz Kiehn
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/// \file
/// \brief   Efficient evaluation of polynomial functions
/// \author  Moritz Kiehn <msmk@cern.ch>
/// \date    2018-02-26

#include <array>
#include <iterator>
#include <type_traits>
#include <utility>

#pragma once

namespace dfe {

/// Evaluate a polynomial of arbitrary order.
///
/// \param x      Where to evaluate the polynomial.
/// \param coeffs ReversibleContainer with n+1 coefficients.
///
/// The coefficients must be given in increasing order. E.g. a second order
/// polynomial has three coefficients c0, c1, c2 that define the function
///
///     f(x) = c0 + c1*x + c2*x^2
///
template<typename T, typename Container>
constexpr T
polynomial_val(const T& x, const Container& coeffs) {
  // Use Horner's method to evaluate polynomial, i.e. expand
  //   f(x) = c0 + c1*x + c2*x^2 + c3*x^3
  // to the following form
  //   f(x) = c0 + x * (c1 + x * (c2 + x * c3))
  // to allow iterative computation with minimal number of operations
  T value = x; // make sure dynamic-sized types, e.g. std::valarray, work
  value = 0;
  for (auto c = std::rbegin(coeffs); c != std::rend(coeffs); ++c) {
    value = *c + x * value;
  }
  return value;
}

/// Evaluate the value and the derivative of a polynomial of arbitrary order.
///
/// \param x      Where to evaluate the derivative.
/// \param coeffs ReversibleContainer with n+1 coefficients.
///
/// The coefficients must be given in increasing order. E.g. a second order
/// polynomial has three coefficients c0, c1, c2 that define the function
///
///     f(x) = c0 + c1*x + c2*x^2
///
/// and the derivative
///
///     df(x)/dx = 0 + c1 + 2*c2*x
///
template<typename T, typename Container>
constexpr std::pair<T, T>
polynomial_valder(const T& x, const Container& coeffs) {
  // Use Horner's method to evaluate polynomial and its derivative at the
  // same time
  T q = x; // make sure dynamic-sized types, e.g. std::valarray, work
  T p = x;
  q = 0;
  p = 0;
  for (auto c = std::rbegin(coeffs); c != std::rend(coeffs); ++c) {
    q = p + x * q;
    p = *c + x * p;
  }
  return {p, q};
}

/// Evaluate the the derivative of a polynomial of arbitrary order.
///
/// \param x      Where to evaluate the derivative.
/// \param coeffs ReversibleContainer with n+1 coefficients.
///
/// The coefficients must be given in increasing order. E.g. a second order
/// polynomial has three coefficients c0, c1, c2 that define the derivative
///
///     df(x)/dx = 0 + c1 + 2*c2*x
///
template<typename T, typename Container>
constexpr T
polynomial_der(const T& x, const Container& coeffs) {
  return polynomial_valder(x, coeffs).second;
}

/// Evaluate a polynomial with an order fixed at compile time.
template<typename T, typename U>
constexpr auto
polynomial_val(const T& x, std::initializer_list<U> coeffs) {
  return polynomial_val<T, std::initializer_list<U>>(x, coeffs);
}

/// Evaluate the derivative of a polynomial with an order fixed at compile time.
template<typename T, typename U>
constexpr auto
polynomial_der(const T& x, std::initializer_list<U> coeffs) {
  return polynomial_der<T, std::initializer_list<U>>(x, coeffs);
}

/// Evaluate the derivative of a polynomial with an order fixed at compile time.
template<typename T, typename U>
constexpr auto
polynomial_valder(const T& x, std::initializer_list<U> coeffs) {
  return polynomial_valder<T, std::initializer_list<U>>(x, coeffs);
}

} // namespace dfe

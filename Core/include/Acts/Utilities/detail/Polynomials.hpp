// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/ArrayHelpers.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <cassert>

namespace Acts::detail {
/// @brief Transforms the coefficients of a n-th degree polynomial
///        into the one corresponding to its D-th derivative
///        The result are (N-D) non-vanishing coefficients
/// @tparam D: Order of the derivative to calculate
/// @tparam N: Order of the original polynomial
/// @param coeffs: Reference to the polynomial's coefficients
template <std::size_t D, std::size_t N>
constexpr std::array<double, N - D> derivativeCoefficients(
    const std::array<double, N>& coeffs) {
  static_assert(N > D, "Coefficients trivially collapse to 0.");

  if constexpr (D == 0) {
    return coeffs;
  } else {
    std::array<double, N - 1> newCoeffs{filledArray<double, N - 1>(0.)};

    for (std::size_t i = 0; i < N - 1; ++i) {
      newCoeffs[i] = (i + 1) * coeffs[i + 1];
    }

    return derivativeCoefficients<D - 1>(newCoeffs);
  }
}

/// @brief Evaluates a polynomial with degree (N-1) at domain value x
/// @param x: Value where the polynomial is to be evaluated
/// @param coeff: Polynomial coefficients, where the i-th entry
///               is associated to the x^{i} monomial
template <std::size_t N>
constexpr double polynomialSum(const double x,
                               const std::array<double, N>& coeffs) {
  double result{0.};
  double y{1.};
  for (unsigned k = 0; k < N; ++k) {
    if (abs(coeffs[k]) > std::numeric_limits<double>::epsilon()) {
      result += coeffs[k] * y;
    }
    y *= x;
  }
  return result;
}

/// @brief Evaluates the D-th derivative of a polynomial with degree (N-1)
///        at domain value x
/// @tparam D: Order of the derivative to calculate
/// @tparam N: Order of the original polynomial
/// @param x: Value where the polynomial is to be evaluated
/// @param coeff: Polynomial coefficients, where the i-th entry
///               is associated to the x^{i} monomial
template <std::size_t D, std::size_t N>
constexpr double derivativeSum(const double x,
                               const std::array<double, N>& coeffs) {
  if constexpr (N > D) {
    return polynomialSum(x, derivativeCoefficients<D, N>(coeffs));
  }
  return 0.;
}

namespace Legendre {
/// @brief Calculates the n-th coefficient of the legendre polynomial series
///        (cf. https://en.wikipedia.org/wiki/Legendre_polynomials)
/// @param l: Order of the legendre polynomial
/// @param k: Coefficient inside the polynomial representation
constexpr double coeff(const unsigned l, const unsigned k) {
  assert(k <= l);
  if ((k % 2) != (l % 2)) {
    return 0.;
  } else if (k > 1) {
    const double a_k = -(1. * (l - k + 2) * (l + k - 1)) /
                       (1. * (k * (k - 1))) * coeff(l, k - 2);
    return a_k;
  }
  unsigned fl = (l - l % 2) / 2;
  unsigned binom = binomial(l, fl) * binomial(2 * l - 2 * fl, l);
  return ((fl % 2) != 0 ? -1. : 1.) * pow(0.5, l) * (1. * binom);
}
/// @brief Assembles the coefficients of the L-th legendre polynomial
///        and returns them via an array
/// @tparam L: Order of the legendre polynomial
template <unsigned L>
constexpr std::array<double, L + 1> coefficients() {
  std::array<double, L + 1> allCoeffs{filledArray<double, L + 1>(0.)};
  for (unsigned k = (L % 2); k <= L; k += 2) {
    allCoeffs[k] = coeff(L, k);
  }
  return allCoeffs;
}
/// @brief Evaluates the Legendre polynomial
/// @param x: Domain value to evaluate
/// @param l: Order of the Legendre polynomial
/// @param d: Order of the derivative
constexpr double evaluate(const double x, const unsigned l, unsigned d = 0u) {
  double sum{0.};
  for (unsigned k = l % 2; k + d <= l; k += 2u) {
    sum += pow(x, k - d) * coeff(l, k) *
           (d > 0u ? factorial(k) / factorial(k - d) : 1u);
  }
  return sum;
}
}  // namespace Legendre

/// @brief Implementation of the chebychev polynomials of the first & second kind
///        (c.f. https://en.wikipedia.org/wiki/Chebyshev_polynomials)
namespace Chebychev {
//// @brief Calculates the k-th coefficient of the Chebychev polynomial of the first kind (T_{n})
/// @param n: Order of the polynomial
/// @param k: Coefficient mapped to the monomial n-2k
constexpr double coeffTn(const unsigned n, const unsigned k) {
  assert(n >= 2 * k);
  if (n == 0) {
    return 1.;
  }
  const double sign = (k % 2 == 1 ? -1. : 1.);
  const double t_k = sign * static_cast<double>(factorial(n - k - 1)) /
                     static_cast<double>(factorial(k) * factorial(n - 2 * k)) *
                     static_cast<double>(n) *
                     pow(2., static_cast<int>(n - 2 * k - 1));
  return t_k;
}
/// @brief Collects the coefficients of the n-th Chebychev polynomial
///        of the first kind and returns them via an array
/// @tparam Tn: Order of the Chebychev 1-st kind polynomial
template <unsigned Tn>
constexpr std::array<double, Tn + 1> coeffientsFirstKind() {
  std::array<double, Tn + 1> coeffs{filledArray<double, Tn + 1>(0.)};
  for (unsigned k = 0; 2 * k <= Tn; ++k) {
    coeffs[Tn - 2 * k] = coeffTn(Tn, k);
  }
  return coeffs;
}
/// @brief Evaluates the chebychev polynomial of the first kind
/// @param x: Domain value to evaluate
/// @param n: Order of the chebychev polynomial
/// @param d: Order of the derivative
constexpr double evalFirstKind(const double x, const unsigned n,
                               const unsigned d = 0u) {
  double result{0.};
  for (unsigned k = 0u; 2u * k + d <= n; ++k) {
    result += coeffTn(n, k) * pow(x, n - 2u * k - d) *
              (d > 0 ? factorial(n - 2u * k) / factorial(n - 2u * k - d) : 1);
  }
  return result;
}
/// @brief Calculates the k-th coefficient of the Chebychev polynomial of the second kind (U_{n})
/// @param n: Order of the polynomial
/// @param k: Coefficient mapped to the monomial n-2k
constexpr double coeffUn(const unsigned n, const unsigned k) {
  assert(n >= 2u * k);
  if (n == 0u) {
    return 1.;
  }
  const double sign = (k % 2 == 1u ? -1. : 1.);
  const double u_k = sign * binomial(n - k, k) * pow(2., n - 2u * k);
  return u_k;
}

/// @brief Collects the coefficients of the n-th Chebychev polynomial
///        of the second kind and returns them via an array
/// @tparam Un: Order of the Chebychev 2-nd kind polynomial
template <unsigned Un>
constexpr std::array<double, Un + 1> coeffientsSecondKind() {
  std::array<double, Un + 1> coeffs{filledArray<double, Un + 1>(0.)};
  for (unsigned k = 0u; 2u * k <= Un; ++k) {
    coeffs[Un - 2u * k] = coeffUn(Un, k);
  }
  return coeffs;
}
/// @brief Evaluates the chebychev polynomial of the second kind
/// @param x: Domain value to evaluate
/// @param n: Order of the chebychev polynomial
/// @param d: Order of the derivative
constexpr double evalSecondKind(const double x, const unsigned n,
                                const unsigned d) {
  double result{0.};
  for (unsigned k = 0u; 2u * k + d <= n; ++k) {
    result += coeffUn(n, k) * pow(x, n - 2u * k - d) *
              (d > 0u ? factorial(n - 2u * k) / factorial(n - 2u * k - d) : 1u);
  }
  return result;
}
}  // namespace Chebychev

/// @brief Helper macros to setup the evaluation of the n-th orthogonal
///        polynomial and of its derivatives. The  coefficients of the
///        n-th polynomial and of the first two derivatives
///        are precalculated at compile time. For higher order derivatives
///        a runtime evaluation function is used
/// @param order: Explicit order of the orthogonal polynomial
/// @param derivative: Requested derivative from the function call
/// @param coeffGenerator: Templated function which returns the coefficients
///                        as an array
#define SETUP_POLY_ORDER(order, derivative, coeffGenerator) \
  case order: {                                             \
    constexpr auto polyCoeffs = coeffGenerator<order>();    \
    switch (derivative) {                                   \
      case 0u:                                              \
        return derivativeSum<0>(x, polyCoeffs);             \
      case 1u:                                              \
        return derivativeSum<1>(x, polyCoeffs);             \
      case 2u:                                              \
        return derivativeSum<2>(x, polyCoeffs);             \
      default:                                              \
        break;                                              \
    }                                                       \
    break;                                                  \
  }
/// @brief Helper macro to write a function to evaluate an orthogonal
///        polynomial. The first 10 polynomial terms together with their
///        first two derivatives are compiled as constexpr expressions
///        for the remainders the fall back runTime evaluation function is used
/// @param polyFuncName: Name of the final function
/// @param coeffGenerator: Name of the templated function over the polynomial order
///                        generating the coefficients of the n-th orthogonal
///                        polynomial
/// @param runTimeEval: Evaluation function with the signature
///                       double runTimeEval(const double x,
///                                          const unsigned order,
///                                          onst unsigned derivative)
///                     which is used as fallback if higher order derivatives or
///                     higher order polynomials are requested from the user.
#define EVALUATE_POLYNOMIAL(polyFuncName, coeffGenerator, runTimeEval) \
  constexpr double polyFuncName(double x, unsigned order,              \
                                unsigned derivative = 0) {             \
    switch (order) {                                                   \
      SETUP_POLY_ORDER(0u, derivative, coeffGenerator)                 \
      SETUP_POLY_ORDER(1u, derivative, coeffGenerator)                 \
      SETUP_POLY_ORDER(2u, derivative, coeffGenerator)                 \
      SETUP_POLY_ORDER(3u, derivative, coeffGenerator)                 \
      SETUP_POLY_ORDER(4u, derivative, coeffGenerator)                 \
      SETUP_POLY_ORDER(5u, derivative, coeffGenerator)                 \
      SETUP_POLY_ORDER(6u, derivative, coeffGenerator)                 \
      SETUP_POLY_ORDER(7u, derivative, coeffGenerator)                 \
      SETUP_POLY_ORDER(8u, derivative, coeffGenerator)                 \
      SETUP_POLY_ORDER(9u, derivative, coeffGenerator)                 \
      SETUP_POLY_ORDER(10u, derivative, coeffGenerator)                \
      default:                                                         \
        break;                                                         \
    }                                                                  \
    return runTimeEval(x, order, derivative);                          \
  }

EVALUATE_POLYNOMIAL(chebychevPolyTn, Chebychev::coeffientsFirstKind,
                    Chebychev::evalFirstKind);
EVALUATE_POLYNOMIAL(chebychevPolyUn, Chebychev::coeffientsSecondKind,
                    Chebychev::evalSecondKind);
EVALUATE_POLYNOMIAL(legendrePoly, Legendre::coefficients, Legendre::evaluate);
#undef EVALUATE_POLYNOMIAL
#undef SETUP_POLY_ORDER

}  // namespace Acts::detail

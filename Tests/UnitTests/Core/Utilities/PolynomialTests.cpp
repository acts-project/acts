// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/detail/Polynomials.hpp"

#include <algorithm>
#include <numeric>
namespace {
constexpr double stepSize = 1.e-4;

constexpr bool withinTolerance(const double value, const double expect) {
  constexpr double tolerance = 1.e-10;
  return std::abs(value - expect) < tolerance;
}
template <typename T, std::size_t D>
std::ostream& operator<<(std::ostream& ostr,

                         const std::array<T, D>& arr

) {
  ostr << "[";
  for (std::size_t t = 0; t < D; ++t) {
    ostr << arr[t];
    if (t + 1 != D) {
      ostr << ", ";
    }
  }
  ostr << "]";
  return ostr;
}

template <std::size_t O, std::size_t D>
bool checkDerivative(const std::array<double, O>& unitArray,
                     const std::array<double, (O - D + 1)>& recentDeriv) {
  const auto nthDerivative = Acts::detail::derivativeCoefficients<D>(unitArray);
  const auto firstDeriv = Acts::detail::derivativeCoefficients<1>(recentDeriv);
  std::cout << D << "-th derivative: " << nthDerivative << std::endl;
  if (nthDerivative.size() != firstDeriv.size()) {
    std::cout << __func__ << "<" << D << ">" << ":" << __LINE__
              << " - nthDerivative: " << nthDerivative.size()
              << " vs. firstDeriv: " << firstDeriv.size() << std::endl;
    return false;
  }
  bool good{true};
  for (std::size_t i = 0; i < nthDerivative.size(); ++i) {
    if (!withinTolerance(firstDeriv[i], nthDerivative[i])) {
      std::cout << __func__ << "<" << D << ">" << ":" << __LINE__
                << " - nthDerivative: " << nthDerivative[i]
                << " vs. firstDeriv: " << firstDeriv[i] << std::endl;
      good = false;
    }
  }
  for (std::size_t i = 0; i < firstDeriv.size(); ++i) {
    if (!withinTolerance(firstDeriv[i], (i + 1) * recentDeriv[i + 1])) {
      std::cout << __func__ << "<" << D << ">" << ":" << __LINE__
                << " - firstDeriv: " << firstDeriv[i]
                << " vs. expect: " << ((i + 1) * recentDeriv[i + 1])
                << std::endl;
      good = false;
    }
  }
  if (!good) {
    return false;
  }
  if constexpr (D + 1 < O) {
    if (!checkDerivative<O, D + 1>(unitArray, firstDeriv)) {
      return false;
    }
  }
  return good;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(PolynomialTests)

BOOST_AUTO_TEST_CASE(DerivativeCoeffs) {
  constexpr std::size_t order = 20;
  std::array<double, order> unitCoeffs{Acts::filledArray<double, order>(1)};
  const bool result = checkDerivative<order, 1>(unitCoeffs, unitCoeffs);
  BOOST_CHECK_EQUAL(result, true);
}
BOOST_AUTO_TEST_CASE(PowerTests) {
  for (unsigned p = 0; p <= 15; ++p) {
    BOOST_CHECK_EQUAL(std::pow(2., p), Acts::pow(2., p));
    BOOST_CHECK_EQUAL(std::pow(0.5, p), Acts::pow(0.5, p));
    for (std::size_t k = 1; k <= 15; ++k) {
      BOOST_CHECK_EQUAL(std::pow(k, p), Acts::pow(k, p));
    }
  }
  for (int p = 0; p <= 15; ++p) {
    BOOST_CHECK_EQUAL(std::pow(2., p), Acts::pow(2., p));
    BOOST_CHECK_EQUAL(std::pow(2., -p), Acts::pow(2., -p));
    BOOST_CHECK_EQUAL(std::pow(0.5, p), Acts::pow(0.5, p));

    BOOST_CHECK_EQUAL(std::pow(0.5, -p), Acts::pow(0.5, -p));
  }
}

BOOST_AUTO_TEST_CASE(SumOfIntegers) {
  std::array<unsigned, 100> numberSeq{Acts::filledArray<unsigned, 100>(1)};
  std::iota(numberSeq.begin(), numberSeq.end(), 1);
  for (unsigned i = 1; i <= numberSeq.size(); ++i) {
    const unsigned sum =
        std::accumulate(numberSeq.begin(), numberSeq.begin() + i, 0);
    BOOST_CHECK_EQUAL(sum, Acts::sumUpToN(i));
  }
}

BOOST_AUTO_TEST_CASE(BinomialTests) {
  BOOST_CHECK_EQUAL(Acts::binomial(1u, 1u), 1u);
  for (unsigned n = 2; n <= 10; ++n) {
    /// Check that the binomial of (n 1 is always n)
    BOOST_CHECK_EQUAL(Acts::binomial(n, 1u), n);
    for (unsigned k = 1; k <= n; ++k) {
      /// Use recursive formula
      ///  n      n -1       n -1
      ///     =          +
      ///  k      k -1        k
      std::cout << "n: " << n << ", k: " << k
                << ", binom(n,k): " << Acts::binomial(n, k)
                << ", binom(n-1, k-1): " << Acts::binomial(n - 1, k - 1)
                << ", binom(n-1,k): " << Acts::binomial(n - 1, k) << std::endl;
      BOOST_CHECK_EQUAL(Acts::binomial(n, k), Acts::binomial(n - 1, k - 1) +
                                                  Acts::binomial(n - 1, k));
      BOOST_CHECK_EQUAL(Acts::binomial(n, k), Acts::binomial(n, n - k));
    }
  }
}
BOOST_AUTO_TEST_CASE(LegendrePolynomials) {
  using namespace Acts::detail::Legendre;
  using namespace Acts::detail;
  std::cout << "Legdnre coefficients L=0: " << coefficients<0>() << std::endl;
  std::cout << "Legdnre coefficients L=1: " << coefficients<1>() << std::endl;
  std::cout << "Legdnre coefficients L=2: " << coefficients<2>() << std::endl;
  std::cout << "Legdnre coefficients L=3: " << coefficients<3>() << std::endl;
  std::cout << "Legdnre coefficients L=4: " << coefficients<4>() << std::endl;
  std::cout << "Legdnre coefficients L=5: " << coefficients<5>() << std::endl;
  std::cout << "Legdnre coefficients L=6: " << coefficients<6>() << std::endl;
  for (unsigned order = 0; order < 10; ++order) {
    const double sign = (order % 2 == 0 ? 1. : -1.);
    BOOST_CHECK_EQUAL(withinTolerance(legendrePoly(1., order), 1.), true);
    BOOST_CHECK_EQUAL(withinTolerance(legendrePoly(-1., order), sign * 1.),
                      true);
    for (double x = -1.; x <= 1.; x += stepSize) {
      const double evalX = legendrePoly(x, order);
      const double evalDx = legendrePoly(x, order, 1);
      const double evalD2x = legendrePoly(x, order, 2);
      /// Check whether the polynomial solves the legendre differental equation
      /// (1-x^{2}) d^P_{l}(x) /d^{x} -2x * d^P_{l}(x) / dx + l*(l+1) P_{l}(x) =
      /// 0;
      const double legendreEq = (1. - Acts::square(x)) * evalD2x -
                                2. * x * evalDx + order * (order + 1) * evalX;
      BOOST_CHECK_EQUAL(withinTolerance(legendreEq, 0.), true);
      BOOST_CHECK_EQUAL(withinTolerance(evalX, sign * legendrePoly(-x, order)),
                        true);
      if (order == 0) {
        continue;
      }
      const double evalX_P1 = legendrePoly(x, order + 1);
      const double evalX_M1 = legendrePoly(x, order - 1);
      /// Recursion formular
      /// (n+1) P_{n+1}(x) = (2n+1)*x*P_{n}(x) - n * P_{n-1}(x)
      BOOST_CHECK_EQUAL(
          withinTolerance((order + 1) * evalX_P1,
                          (2 * order + 1) * x * evalX - order * evalX_M1),
          true);
    }
  }
}

BOOST_AUTO_TEST_CASE(ChebychevPolynomials) {
  using namespace Acts::detail::Chebychev;
  using namespace Acts::detail;

  /// Properties of the Chebychev polynomials
  ///  T_{n}(x) = (-1)^{n} T_n(-x)
  ///  U_{n}(x) = (-1)^{n} U_n(-x)
  ///  T_{n}(1) = 1
  ///  U_{n}(1) = n+1
  ///  T_{n+1}(x) = 2*x*T_{n}(x) - T_{n-1}(x)
  ///             = x*T_{n}(x) - (1-x^{2}) U_{n-1}(x)
  ///  U_{n+1}(x) = 2*x*U_{n}(x) - U_{n-1}(x)
  ///             = x*U_{n}(x) + T_{n+1}(x)
  for (unsigned order = 0; order < 10; ++order) {
    const double sign = (order % 2 == 0 ? 1. : -1.);
    const double T_n1 = chebychevPolyTn(1., order);
    const double U_n1 = chebychevPolyUn(1., order);

    std::cout << "Order: " << order << " T(1)=" << T_n1 << ", U(1)=" << U_n1
              << std::endl;
    BOOST_CHECK_EQUAL(withinTolerance(T_n1, 1.), true);
    BOOST_CHECK_EQUAL(withinTolerance(chebychevPolyTn(-1., order), sign), true);

    BOOST_CHECK_EQUAL(withinTolerance(U_n1, order + 1), true);
    BOOST_CHECK_EQUAL(withinTolerance(U_n1, sign * chebychevPolyUn(-1., order)),
                      true);
    if (order == 0) {
      continue;
    }
    for (double x = -1.; x <= 1.; x += stepSize) {
      const double U_np1 = chebychevPolyUn(x, order + 1);
      const double U_n = chebychevPolyUn(x, order);
      const double U_nm1 = chebychevPolyUn(x, order - 1);

      const double T_np1 = chebychevPolyTn(x, order + 1);
      const double T_n = chebychevPolyTn(x, order);
      const double T_nm1 = chebychevPolyTn(x, order - 1);

      BOOST_TEST_MESSAGE(
          "Order: " << order << ", x=" << x << ", U_{n+1}(x) = " << U_np1
                    << ", 2*x*U_{n}(x) - U_{n-1}=" << (2. * x * U_n - U_nm1)
                    << ", x*U_{n}(x) + T_{n+1}(x)=" << (x * U_n + T_np1));

      BOOST_CHECK_EQUAL(withinTolerance(U_np1, 2. * x * U_n - U_nm1), true);
      BOOST_CHECK_EQUAL(withinTolerance(U_np1, x * U_n + T_np1), true);

      BOOST_CHECK_EQUAL(withinTolerance(U_n, sign * chebychevPolyUn(-x, order)),
                        true);

      BOOST_TEST_MESSAGE("x=" << x << ", T_{n+1}(x) = " << T_np1
                              << ", 2*x*T_{n}(x) - T_{n-1}="
                              << (2. * x * T_n - T_nm1));
      BOOST_CHECK_EQUAL(withinTolerance(T_np1, 2. * x * T_n - T_nm1), true);
      BOOST_CHECK_EQUAL(withinTolerance(T_np1, x * T_n - (1 - x * x) * U_nm1),
                        true);

      BOOST_CHECK_EQUAL(withinTolerance(T_n, sign * chebychevPolyTn(-x, order)),
                        true);

      /// Check derivative
      BOOST_CHECK_EQUAL(withinTolerance(chebychevPolyTn(x, order, 1),
                                        Chebychev::evalFirstKind(x, order, 1)),
                        true);
      BOOST_CHECK_EQUAL(withinTolerance(chebychevPolyUn(x, order, 1),
                                        Chebychev::evalSecondKind(x, order, 1)),
                        true);
      BOOST_CHECK_EQUAL(withinTolerance(chebychevPolyTn(x, order, 2),
                                        Chebychev::evalFirstKind(x, order, 2)),
                        true);
      BOOST_CHECK_EQUAL(withinTolerance(chebychevPolyUn(x, order, 2),
                                        Chebychev::evalSecondKind(x, order, 2)),
                        true);
    }
  }
}
BOOST_AUTO_TEST_SUITE_END()

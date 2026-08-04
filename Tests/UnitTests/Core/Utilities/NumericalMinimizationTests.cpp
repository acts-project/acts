// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/detail/NumericalMinimization.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <limits>
#include <string>

using namespace Acts;
using namespace Acts::detail;

namespace {

/// An objective that ScalarObjective must accept
struct GoodObjective {
  double operator()(const Vector<2>& p) const { return p.squaredNorm(); }
};

/// An objective ScalarObjective must reject: not callable with a Vector<2>
struct BadObjective {
  double operator()(const std::string& /*p*/) const { return 0.0; }
};

static_assert(ScalarObjective<GoodObjective, 2>);
static_assert(!ScalarObjective<BadObjective, 2>);

}  // namespace

BOOST_AUTO_TEST_SUITE(NumericalMinimizationSuite)

BOOST_AUTO_TEST_CASE(NelderMead_Quadratic1D) {
  const Vector<1> minimum{2.5};
  const auto objective = [&minimum](const Vector<1>& p) {
    return (p - minimum).squaredNorm();
  };

  const auto result = nelderMead<1>(objective, Vector<1>{0.0}, Vector<1>{1.0});
  BOOST_REQUIRE(result.ok());
  CHECK_CLOSE_ABS((*result)(0), minimum(0), 1e-6);
}

BOOST_AUTO_TEST_CASE(NelderMead_Quadratic2D) {
  const Vector<2> minimum{1.2, -3.4};
  const auto objective = [&minimum](const Vector<2>& p) {
    return (p - minimum).squaredNorm();
  };

  const auto result =
      nelderMead<2>(objective, Vector<2>{0.0, 0.0}, Vector<2>{1.0, 1.0});
  BOOST_REQUIRE(result.ok());
  CHECK_CLOSE_ABS((*result)(0), minimum(0), 1e-6);
  CHECK_CLOSE_ABS((*result)(1), minimum(1), 1e-6);
}

BOOST_AUTO_TEST_CASE(NelderMead_Quadratic3D) {
  const Vector<3> minimum{-1.0, 2.0, 0.5};
  const auto objective = [&minimum](const Vector<3>& p) {
    return (p - minimum).squaredNorm();
  };

  const auto result = nelderMead<3>(objective, Vector<3>{0.0, 0.0, 0.0},
                                    Vector<3>{1.0, 1.0, 1.0});
  BOOST_REQUIRE(result.ok());
  CHECK_CLOSE_ABS((*result)(0), minimum(0), 1e-6);
  CHECK_CLOSE_ABS((*result)(1), minimum(1), 1e-6);
  CHECK_CLOSE_ABS((*result)(2), minimum(2), 1e-6);
}

BOOST_AUTO_TEST_CASE(NelderMead_Rosenbrock2D) {
  // Standard Rosenbrock banana, minimum at (1, 1)
  const auto objective = [](const Vector<2>& p) {
    const double a = 1 - p(0);
    const double b = p(1) - p(0) * p(0);
    return a * a + 100 * b * b;
  };

  const auto result =
      nelderMead<2>(objective, Vector<2>{-1.0, 1.0}, Vector<2>{0.5, 0.5});
  BOOST_REQUIRE(result.ok());
  CHECK_CLOSE_ABS((*result)(0), 1.0, 1e-3);
  CHECK_CLOSE_ABS((*result)(1), 1.0, 1e-3);
}

BOOST_AUTO_TEST_CASE(NelderMead_AllInfinite_Fails) {
  const auto objective = [](const Vector<2>& /*p*/) {
    return std::numeric_limits<double>::infinity();
  };

  const auto result =
      nelderMead<2>(objective, Vector<2>{0.0, 0.0}, Vector<2>{1.0, 1.0});
  BOOST_CHECK(!result.ok());
}

BOOST_AUTO_TEST_CASE(NumericalCovariance_Quadratic2D) {
  // f(x, y) = x^2 / (2 a^2) + y^2 / (2 b^2) has Hessian diag(1/a^2, 1/b^2),
  // so the inverse (the covariance) is diag(a^2, b^2)
  const double a = 2.0;
  const double b = 0.5;
  const auto objective = [a, b](const Vector<2>& p) {
    return 0.5 * p(0) * p(0) / (a * a) + 0.5 * p(1) * p(1) / (b * b);
  };

  const auto covariance = numericalCovariance<2>(objective, Vector<2>{0.0, 0.0},
                                                 Vector<2>{0.1 * a, 0.1 * b});
  BOOST_REQUIRE(covariance.ok());
  CHECK_CLOSE_REL((*covariance)(0, 0), a * a, 1e-6);
  CHECK_CLOSE_REL((*covariance)(1, 1), b * b, 1e-6);
  CHECK_CLOSE_ABS((*covariance)(0, 1), 0.0, 1e-9);
  CHECK_CLOSE_ABS((*covariance)(1, 0), 0.0, 1e-9);
}

BOOST_AUTO_TEST_CASE(NumericalCovariance_NotAMinimum_Fails) {
  // Saddle point: Hessian is diag(2, -2), not positive definite
  const auto objective = [](const Vector<2>& p) {
    return p(0) * p(0) - p(1) * p(1);
  };

  const auto covariance = numericalCovariance<2>(objective, Vector<2>{0.0, 0.0},
                                                 Vector<2>{0.1, 0.1});
  BOOST_CHECK(!covariance.ok());
}

BOOST_AUTO_TEST_CASE(NumericalCovariance_NonPositiveStep_Fails) {
  const auto objective = [](const Vector<2>& p) { return p.squaredNorm(); };

  const auto covariance = numericalCovariance<2>(objective, Vector<2>{0.0, 0.0},
                                                 Vector<2>{0.1, 0.0});
  BOOST_CHECK(!covariance.ok());
}

BOOST_AUTO_TEST_SUITE_END()

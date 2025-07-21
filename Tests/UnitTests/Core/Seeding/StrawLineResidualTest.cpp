// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/detail/StrawLineFitAuxiliaries.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
using namespace Acts;

using namespace Acts::detail;
using namespace Acts::UnitLiterals;
using Line_t = StrawLineFitAuxiliaries::Line_t;
using Vector = Line_t::Vector;
using Config_t = StrawLineFitAuxiliaries::Config;
using Pars_t = Line_t::ParamVector;
using ParIdx = Line_t::ParIndices;

namespace {
std::string parName(const unsigned idx) {
  switch (idx) {
    using enum ParIdx;
    case x0:
      return "x0";
    case y0:
      return "y0";
    case theta:
      return "theta";
    case phi:
      return "phi";
    default:
      break;
  }
  return "unknown";
}

}  // namespace
namespace Acts::Test {
class TestSpacePoint {
 public:
  TestSpacePoint(
      const Vector3& pos, const Vector3& wireDir, const double radius,
      const ActsSquareMatrix<3>& cov = ActsSquareMatrix<3>::Identity())
      : m_pos{pos}, m_dir{wireDir}, m_radius{radius}, m_cov{cov} {}

  TestSpacePoint(
      const Vector3& pos, const Vector3& stripDir, const Vector3& stripNorm,
      const ActsSquareMatrix<3>& cov = ActsSquareMatrix<3>::Identity())
      : m_pos{pos},
        m_dir{stripDir},
        m_toNext{stripNorm},
        m_cov{cov},
        m_isStraw{false} {}

  const Vector3& localPosition() const { return m_pos; }
  const Vector3& sensorDirection() const { return m_dir; }
  const Vector3& sensorNormal() const { return m_toNext; }
  const Vector3& planeNormal() const { return m_planeNorm; }

  bool isStraw() const { return m_isStraw; }

  double driftRadius() const { return m_radius; }
  double time() const { return m_time; }
  const ActsSquareMatrix<3>& covariance() const { return m_cov; }

 private:
  Vector3 m_pos{Vector3::Zero()};
  Vector3 m_dir{Vector3::Zero()};
  Vector3 m_toNext{Vector3::Zero()};
  Vector3 m_planeNorm{m_dir.cross(m_toNext)};
  double m_radius{0.};
  double m_time{0.};
  ActsSquareMatrix<3> m_cov{ActsSquareMatrix<3>::Zero()};
  bool m_isStraw{true};
};

BOOST_AUTO_TEST_SUITE(StrawLineSeederTest)

void testResidual(const Pars_t& linePars, const TestSpacePoint& testPoint) {
  Config_t resCfg{};
  resCfg.useHessian = true;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::y0, ParIdx::theta, ParIdx::phi};
  std::cout << "Test residual w.r.t. sp with position: "
            << toString(testPoint.localPosition())
            << ", orientation: " << toString(testPoint.sensorDirection())
            << ", r: " << testPoint.driftRadius() << std::endl;
  Line_t line{};
  line.updateParameters(linePars);

  std::cout << "WireResidualTest - Test segment: " << toString(line.position())
            << ", " << toString(line.direction()) << std::endl;

  StrawLineFitAuxiliaries resCalc{resCfg};
  if (testPoint.isStraw()) {
    resCalc.updateStrawResidual(line, testPoint);
  } else {
    resCalc.updateStripResidual(line, testPoint);
  }
  using namespace Acts::detail::LineHelper;
  const double lineDist =
      signedDistance(testPoint.localPosition(), testPoint.sensorDirection(),
                     line.position(), line.direction());

  std::cout << "Residual: " << toString(resCalc.residual())
            << ", line distance: " << lineDist << std::endl;
  ///
  BOOST_CHECK_CLOSE(resCalc.residual().norm(),
                    lineDist - testPoint.driftRadius(), 1.e-12);
  constexpr double h = 5.e-8;
  constexpr double tolerance = 1.e-5;
  for (auto par : resCfg.parsToUse) {
    Pars_t lineParsUp{linePars}, lineParsDn{linePars};
    lineParsUp[par] += h;
    lineParsDn[par] -= h;
    Line_t lineUp{}, lineDn{};
    lineUp.updateParameters(lineParsUp);
    lineDn.updateParameters(lineParsDn);

    StrawLineFitAuxiliaries resCalcUp{resCfg}, resCalcDn{resCfg};
    if (testPoint.isStraw()) {
      resCalcUp.updateStrawResidual(lineUp, testPoint);
      resCalcDn.updateStrawResidual(lineDn, testPoint);
    } else {
      resCalcUp.updateStripResidual(lineUp, testPoint);
      resCalcDn.updateStripResidual(lineDn, testPoint);
    }
    const Vector numDeriv =
        (resCalcUp.residual() - resCalcDn.residual()) / (2. * h);

    std::cout << "Derivative test: " << parName(par)
              << ", derivative: " << toString(resCalc.gradient(par))
              << ",  numerical: " << toString(numDeriv) << std::endl;
    BOOST_CHECK_LE((numDeriv - resCalc.gradient(par)).norm() /
                       std::max(numDeriv.norm(), 1.),
                   tolerance);
    /// Next attempt to calculate the second derivative
    for (unsigned par1 = 0; par1 <= par; ++par1) {
      Pars_t lineParsUp1{linePars}, lineParsDn1{linePars};
      lineParsUp1[par1] += h;
      lineParsDn1[par1] -= h;

      lineUp.updateParameters(lineParsUp1);
      lineDn.updateParameters(lineParsDn1);
      if (testPoint.isStraw()) {
        resCalcUp.updateStrawResidual(lineUp, testPoint);
        resCalcDn.updateStrawResidual(lineDn, testPoint);
      } else {
        resCalcUp.updateStripResidual(lineUp, testPoint);
        resCalcDn.updateStripResidual(lineDn, testPoint);
      }
      const Vector numDeriv1{
          (resCalcUp.gradient(par) - resCalcDn.gradient(par)) / (2. * h)};
      const Vector& analyticDeriv = resCalc.hessian(par, par1);
      std::cout << "Second deriv (" << parName(par) << ", " << parName(par1)
                << ") -- numerical: " << toString(numDeriv1)
                << ", analytic: " << toString(analyticDeriv) << std::endl;
      BOOST_CHECK_LE((numDeriv1 - analyticDeriv).norm(), tolerance);
    }
  }
}

BOOST_AUTO_TEST_CASE(WireResidualTest) {
  /// Set the line to be 45 degrees
  using Pars_t = Line_t::ParamVector;
  using ParIdx = Line_t::ParIndices;

  Pars_t linePars{};
  linePars[ParIdx::phi] = 90 * 1_degree;
  linePars[ParIdx::theta] = 45 * 1_degree;

  const Vector wireDir = Vector::UnitX();
  /// Generate the first test measurement
  const TestSpacePoint firstTestPoint{Vector{100. * 1_cm, 50. * 1_cm, 0.},
                                      wireDir, 10. * 1_cm};
  testResidual(linePars, firstTestPoint);
}

BOOST_AUTO_TEST_CASE(StripResidual) {
  const Vector dir1{makeDirectionFromPhiTheta(60 * 1_degree, 90 * 1_degree)};
  const Vector dir2{makeDirectionFromPhiTheta(30 * 1_degree, 90 * 1_degree)};

  const Vector pos{474. * 1_cm, -251 * 1_cm, 0.};
  const double dirDots = dir1.dot(dir2);

  ///
  ///       A = lambda * b_{1} + kappa * b_{2}
  ///
  ///                 <A, b_{2}> - <b_{1}, b_{2}> * <A, b_{1}>
  ///   --> kappa  = -----------------------------------------
  ///                         1 -  <b_{1}, b_{2}>^{2}
  ///
  ///                 <A, b_{1}> - <b_{1}, b_{2}> * <A, b_{2}>
  ///   --> lambda = -------------------------------------------
  ///                         1 -  <b_{1}, b_{2}>^{2}

  const double kappa = pos.dot(dir2 - dirDots * dir1) / (1 - square(dirDots));
  const double lambda = pos.dot(dir1 - dirDots * dir2) / (1 - square(dirDots));

  const Vector expand = lambda * dir1 + kappa * dir2;
  std::cout << "Original: " << toString(pos)
            << ", expansion: " << toString(expand) << std::endl;
}

}  // namespace Acts::Test
BOOST_AUTO_TEST_SUITE_END()

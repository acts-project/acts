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
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
using namespace Acts;

using namespace Acts::detail;
using namespace Acts::UnitLiterals;
using namespace Acts::PlanarHelper;
using Line_t = StrawLineFitAuxiliaries::Line_t;
using Vector = Line_t::Vector;
using Config_t = StrawLineFitAuxiliaries::Config;
using Pars_t = Line_t::ParamVector;
using ParIdx = Line_t::ParIndices;

namespace {
std::string parName(const std::size_t idx) {
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
  /// @brief Mask to express which track directions are measured
  ///        by the space point
  enum ProjectorMask : std::size_t {
    bendingDir = 1 << 0,
    nonBendingDir = 1 << 1,
    timeOfArrival = 1 << 2,
    bothDirections = bendingDir | nonBendingDir
  };

  TestSpacePoint(
      const Vector3& pos, const Vector3& wireDir, const double radius,
      bool measNonPrec = false,
      const ActsSquareMatrix<3>& cov = ActsSquareMatrix<3>::Identity())
      : m_pos{pos},
        m_dir{wireDir},
        m_radius{radius},
        m_dirMask{measNonPrec ? bothDirections : bendingDir},
        m_cov{cov} {}

  TestSpacePoint(
      const Vector3& pos, const Vector3& stripDir, const Vector3& stripNorm,
      ProjectorMask mask,
      const ActsSquareMatrix<3>& cov = ActsSquareMatrix<3>::Identity())
      : m_pos{pos},
        m_dir{stripDir},
        m_toNext{stripNorm},
        m_dirMask{mask},
        m_cov{cov},
        m_isStraw{false} {}

  const Vector3& localPosition() const { return m_pos; }
  const Vector3& sensorDirection() const { return m_dir; }
  const Vector3& sensorNormal() const { return m_toNext; }
  const Vector3& planeNormal() const { return m_planeNorm; }

  bool isStraw() const { return m_isStraw; }
  bool hasTime() const { return m_dirMask & (ProjectorMask::timeOfArrival); }
  bool measPrecCoord() const { return m_dirMask & (ProjectorMask::bendingDir); }
  bool measNonPrecCoord() const {
    return m_dirMask & (ProjectorMask::nonBendingDir);
  }

  double driftRadius() const { return m_radius; }
  double time() const { return m_time; }
  const ActsSquareMatrix<3>& covariance() const { return m_cov; }

 private:
  Vector3 m_pos{Vector3::Zero()};
  Vector3 m_dir{Vector3::Zero()};
  Vector3 m_toNext{Vector3::Zero()};
  Vector3 m_planeNorm{m_dir.cross(m_toNext).normalized()};
  double m_radius{0.};
  double m_time{0.};
  ProjectorMask m_dirMask{ProjectorMask::bothDirections};
  ActsSquareMatrix<3> m_cov{ActsSquareMatrix<3>::Zero()};
  bool m_isStraw{true};
};

BOOST_AUTO_TEST_SUITE(StrawLineSeederTest)

void testResidual(const Pars_t& linePars, const TestSpacePoint& testPoint) {
  using namespace Acts::detail::LineHelper;
  using ResidualIdx = StrawLineFitAuxiliaries::ResidualIdx;
  Config_t resCfg{};
  resCfg.useHessian = true;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::y0, ParIdx::phi, ParIdx::theta};
  Line_t line{};
  line.updateParameters(linePars);

  std::cout << "\n\n\nResidual test - Test line: " << toString(line.position())
            << ", " << toString(line.direction()) << std::endl;

  StrawLineFitAuxiliaries resCalc{
      resCfg, Acts::getDefaultLogger("testRes", Logging::Level::VERBOSE)};
  resCalc.updateSpatialResidual(line, testPoint);
  if (testPoint.isStraw()) {
    const double lineDist =
        signedDistance(testPoint.localPosition(), testPoint.sensorDirection(),
                       line.position(), line.direction());
    std::cout << "Test residual w.r.t. sp with position: "
              << toString(testPoint.localPosition())
              << ", orientation: " << toString(testPoint.sensorDirection())
              << ", r: " << testPoint.driftRadius() << std::endl;

    std::cout << "Residual: " << toString(resCalc.residual())
              << ", line distance: " << lineDist << std::endl;
    ///
    BOOST_CHECK_CLOSE(resCalc.residual()[ResidualIdx::bending],
                      lineDist - testPoint.driftRadius(), 1.e-12);
  } else {
    std::cout << "Test residual w.r.t strip space point -- position: "
              << toString(testPoint.localPosition())
              << ", normal: " << toString(testPoint.planeNormal())
              << ", strip: " << toString(testPoint.sensorDirection())
              << ", to-next: " << toString(testPoint.sensorNormal())
              << std::endl;
    const Vector3& n{testPoint.planeNormal()};
    const Vector3& b1 = testPoint.measPrecCoord() ? testPoint.sensorNormal()
                                                  : testPoint.sensorDirection();
    const Vector3& b2 = testPoint.measPrecCoord() ? testPoint.sensorDirection()
                                                  : testPoint.sensorNormal();
    const Vector3 planeIsect = intersectPlane(line.position(), line.direction(),
                                              n, testPoint.localPosition())
                                   .position();
    const Vector3 delta = (planeIsect - testPoint.localPosition());

    std::cout << "Residual: " << toString(resCalc.residual())
              << ", delta: " << toString(delta) << std::endl;
    const Vector3 res = delta - delta.dot(n) * n;
    BOOST_CHECK_LE((res - resCalc.residual()[ResidualIdx::bending] * b1 -
                    resCalc.residual()[ResidualIdx::nonBending] * b2)
                       .norm(),
                   1.e-10);
  }

  constexpr double h = 5.e-9;
  constexpr double tolerance = 1.e-3;
  for (auto par : resCfg.parsToUse) {
    Pars_t lineParsUp{linePars}, lineParsDn{linePars};
    lineParsUp[par] += h;
    lineParsDn[par] -= h;
    Line_t lineUp{}, lineDn{};
    lineUp.updateParameters(lineParsUp);
    lineDn.updateParameters(lineParsDn);

    StrawLineFitAuxiliaries resCalcUp{
        resCfg, Acts::getDefaultLogger("testResUp", Logging::Level::VERBOSE)};
    StrawLineFitAuxiliaries resCalcDn{
        resCfg, Acts::getDefaultLogger("testResDn", Logging::Level::VERBOSE)};

    resCalcUp.updateSpatialResidual(lineUp, testPoint);
    resCalcDn.updateSpatialResidual(lineDn, testPoint);

    const Vector numDeriv =
        (resCalcUp.residual() - resCalcDn.residual()) / (2. * h);

    std::cout << "Derivative test: " << parName(par)
              << ", derivative: " << toString(resCalc.gradient(par))
              << ",  numerical: " << toString(numDeriv) << std::endl;
    BOOST_CHECK_LE((numDeriv - resCalc.gradient(par)).norm() /
                       std::max(numDeriv.norm(), 1.),
                   tolerance);
    /// Next attempt to calculate the second derivative
    for (auto par1 : resCfg.parsToUse) {
      lineParsUp = linePars;
      lineParsDn = linePars;
      lineParsUp[par1] += h;
      lineParsDn[par1] -= h;

      lineUp.updateParameters(lineParsUp);
      lineDn.updateParameters(lineParsDn);

      resCalcUp.updateSpatialResidual(lineUp, testPoint);
      resCalcDn.updateSpatialResidual(lineDn, testPoint);

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
  linePars[ParIdx::phi] = 90. * 1_degree;
  linePars[ParIdx::theta] = 45. * 1_degree;

  const Vector wireDir = Vector::UnitX();
  // /// Generate the first test measurement
  testResidual(linePars,
               TestSpacePoint{Vector{100. * 1_cm, 50. * 1_cm, 30. * 1_cm},
                              wireDir, 10. * 1_cm});
  testResidual(
      linePars,
      TestSpacePoint{Vector{100. * 1_cm, 50. * 1_cm, 30. * 1_cm},
                     makeDirectionFromPhiTheta(10. * 1_degree, 30. * 1_degree),
                     10. * 1_cm});

  testResidual(linePars,
               TestSpacePoint{Vector{100. * 1_cm, 50. * 1_cm, 30. * 1_cm},
                              wireDir, 10. * 1_cm, true});

  linePars[ParIdx::phi] = 30. * 1_degree;
  linePars[ParIdx::theta] = 60. * 1_degree;

  testResidual(linePars,
               TestSpacePoint{Vector{100. * 1_cm, 50. * 1_cm, 30. * 1_cm},
                              wireDir, 10. * 1_cm});

  testResidual(linePars,
               TestSpacePoint{Vector{100. * 1_cm, 50. * 1_cm, 30. * 1_cm},
                              wireDir, 10. * 1_cm, true});

  linePars[ParIdx::phi] = 60. * 1_degree;
  linePars[ParIdx::theta] = 30. * 1_degree;

  testResidual(linePars,
               TestSpacePoint{Vector{100. * 1_cm, 50. * 1_cm, 30. * 1_cm},
                              wireDir, 10. * 1_cm});
  testResidual(linePars,
               TestSpacePoint{Vector{100. * 1_cm, 50. * 1_cm, 30. * 1_cm},
                              wireDir, 10. * 1_cm, true});
}

BOOST_AUTO_TEST_CASE(StripResidual) {
  Pars_t linePars{};
  linePars[ParIdx::phi] = 60. * 1_degree;
  linePars[ParIdx::theta] = 45 * 1_degree;
  return;

  auto placeSp = [&linePars](const double dx, const double dy, const double z) {
    Line_t line{};
    line.updateParameters(linePars);
    const auto pIsect =
        intersectPlane(line.position(), line.direction(), Vector::UnitZ(), z);
    const Acts::Vector3 returnMe =
        (pIsect.position() + Acts::Vector3{dx, dy, 0.});
    return returnMe;
  };

  testResidual(
      linePars,
      TestSpacePoint{placeSp(50. * 1_cm, 50. * 1_cm, 100. * 1_cm),
                     makeDirectionFromPhiTheta(90. * 1_degree, 90 * 1_degree),
                     makeDirectionFromPhiTheta(0. * 1_degree, 90 * 1_degree),
                     TestSpacePoint::bendingDir});

  // testResidual(linePars,
  //                TestSpacePoint{
  //                    pos,
  //                    makeDirectionFromPhiTheta(90. * 1_degree, 90 *
  //                    1_degree), makeDirectionFromPhiTheta(0. * 1_degree, 90*
  //                    1_degree), TestSpacePoint::nonBendingDir});

  // testResidual(linePars,
  //                TestSpacePoint{
  //                    pos,
  //                    makeDirectionFromPhiTheta(0. * 1_degree, 90 * 1_degree),
  //                    makeDirectionFromPhiTheta(90. * 1_degree, 90* 1_degree),
  //                    TestSpacePoint::bothDirections});
  // testResidual(
  //     linePars,
  //     TestSpacePoint{placeSp(50.*1_cm, 50.*1_cm, 100.*1_cm),
  //                    makeDirectionFromPhiTheta(90. * 1_degree, 90 *
  //                    1_degree), makeDirectionFromPhiTheta(0. * 1_degree, 90 *
  //                    1_degree), TestSpacePoint::bendingDir});
}

}  // namespace Acts::Test
BOOST_AUTO_TEST_SUITE_END()

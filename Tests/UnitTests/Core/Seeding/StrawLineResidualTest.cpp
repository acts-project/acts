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
namespace Acts::Test {
class TestSpacePoint {
 public:
  TestSpacePoint(
      const Vector3& pos, const Vector3& dir, const double radius,
      const double time = 0.,
      const ActsSquareMatrix<3>& cov = ActsSquareMatrix<3>::Identity())
      : m_pos{pos}, m_dir{dir}, m_radius{radius}, m_time{time}, m_cov{cov} {}

  const Vector3& localPosition() const { return m_pos; }
  const Vector3& sensorDirection() const { return m_dir; }

  const Vector3& stripPlaneNormal() const { return m_dir; }

  double driftRadius() const { return m_radius; }
  double time() const { return m_time; }
  const ActsSquareMatrix<3>& covariance() const { return m_cov; }

 private:
  Vector3 m_pos{Vector3::Zero()};
  Vector3 m_dir{Vector3::Zero()};
  double m_radius{0.};
  double m_time{0.};
  ActsSquareMatrix<3> m_cov{ActsSquareMatrix<3>::Zero()};
};

Vector3 projectOntoPlane(const Vector3& lineDir, const Vector3& wire) {
  return (lineDir - lineDir.dot(wire) * wire).normalized();
}

BOOST_AUTO_TEST_SUITE(StrawLineSeederTest)

BOOST_AUTO_TEST_CASE(WireResidualTest) {
  using namespace Acts::detail;
  using namespace Acts::UnitLiterals;
  using Line_t = StrawLineFitAuxiliaries::Line_t;
  using Vector = Line_t::Vector;
  using Config_t = StrawLineFitAuxiliaries::Config;
  using Pars_t = Line_t::ParamVector;
  using ParIdx = Line_t::ParIndices;

  /// Set the line to be 45 degrees
  Line_t line{};

  using Pars_t = Line_t::ParamVector;
  using ParIdx = Line_t::ParIndices;
  Pars_t linePars{};
  linePars[ParIdx::phi] = 90 * 1_degree;
  linePars[ParIdx::theta] = 45 * 1_degree;
  line.updateParameters(linePars);

  std::cout << "WireResidualTest - Test segment: " << toString(line.position())
            << ", " << toString(line.direction()) << std::endl;
  const Vector wireDir = Vector::UnitX();
  const Vector conDir = line.direction().cross(wireDir).normalized();
  /// Generate the first test measurement
  constexpr double measDist = 100 * 1_cm;
  constexpr double lineDist = 50. * 1_cm;
  constexpr double driftRadius = 10. * 1_cm;

  const TestSpacePoint firstTestPoint{line.point(measDist) + lineDist * conDir,
                                      wireDir, driftRadius};
  std::cout << "Test point 1:  position: "
            << toString(firstTestPoint.localPosition())
            << ", orientation: " << toString(firstTestPoint.sensorDirection())
            << ", r: " << firstTestPoint.driftRadius() << std::endl;

  Config_t resCfg{};

  StrawLineFitAuxiliaries resCalc{resCfg};
  resCalc.updateStrawResidual(line, firstTestPoint);
  using namespace Acts::detail::LineHelper;
  std::cout << "Residual: " << toString(resCalc.residual())
            << ", line distance: "
            << signedDistance(firstTestPoint.localPosition(),
                              firstTestPoint.sensorDirection(), line.position(),
                              line.direction())
            << std::endl;
  ///
  BOOST_CHECK_CLOSE(resCalc.residual().norm(), lineDist - driftRadius,
                    StrawLineFitAuxiliaries::s_tolerance);
  constexpr double h = 1.e-8;
  constexpr double tolerance = 1.e-7;
  for (auto par : {/*ParIdx::x0, ParIdx::y0,*/ ParIdx::theta, ParIdx::phi}) {
    Pars_t lineParsUp{linePars}, lineParsDn{linePars};
    lineParsUp[par] += h;
    lineParsDn[par] -= h;
    Line_t lineUp{}, lineDn{};
    lineUp.updateParameters(lineParsUp);
    lineDn.updateParameters(lineParsDn);

    StrawLineFitAuxiliaries resCalcUp{resCfg}, resCalcDn{resCfg};
    resCalcUp.updateStrawResidual(lineUp, firstTestPoint);
    resCalcDn.updateStrawResidual(lineDn, firstTestPoint);

    const Vector numDeriv =
        (resCalcUp.residual() - resCalcDn.residual()) / (2. * h);
    const Vector3 projectDeriv =
        (projectOntoPlane(lineUp.direction(), wireDir) -
         projectOntoPlane(lineDn.direction(), wireDir)) /
        (2. * h);
    std::cout << "Derivative test: " << par
              << ", derivative: " << toString(resCalc.gradient(par))
              << ",  numerical: " << toString(numDeriv)
              << ", projection: " << toString(projectDeriv) << std::endl;
    BOOST_CHECK_LE((numDeriv - resCalc.gradient(par)).norm() /
                       std::max(numDeriv.norm(), 1.),
                   tolerance);
    /// Next attempt to calculate the second derivative
  }
}

}  // namespace Acts::Test
BOOST_AUTO_TEST_SUITE_END()

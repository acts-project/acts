// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/CombinatorialSeedSolver.hpp"
#include "Acts/Seeding/detail/StrawLineFitAuxiliaries.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <random>

using namespace Acts;
using namespace Acts::detail;
using namespace Acts::Experimental::detail;
using namespace Acts::UnitLiterals;
using namespace Acts::PlanarHelper;
using RandomEngine = std::mt19937;

using Line_t = StrawLineFitAuxiliaries::Line_t;
using Config_t = StrawLineFitAuxiliaries::Config;
using ParIdx = StrawLineFitAuxiliaries::FitParIndex;

using Vector = Line_t::Vector;
using Pars_t = Line_t::ParamVector;

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
      const std::array<double, 3>& cov = Acts::filledArray<double, 3>(0.))
      : m_pos{pos},
        m_dir{wireDir},
        m_radius{radius},
        m_dirMask{measNonPrec ? bothDirections : bendingDir},
        m_cov{cov} {}

  TestSpacePoint(
      const Vector3& pos, const Vector3& stripDir, const Vector3& stripNorm,
      ProjectorMask mask,
      const std::array<double, 3>& cov = Acts::filledArray<double, 3>(0.))
      : m_pos{pos},
        m_dir{stripDir},
        m_toNext{stripNorm},
        m_dirMask{mask},
        m_cov{cov},
        m_isStraw{false} {}

  const Vector3& localPosition() const { return m_pos; }
  const Vector3& sensorDirection() const { return m_dir; }
  const Vector3& toNextSensor() const { return m_toNext; }
  const Vector3& planeNormal() const { return m_planeNorm; }

  bool isStraw() const { return m_isStraw; }
  bool hasTime() const {
    return (m_dirMask & (ProjectorMask::timeOfArrival)) != 0u;
  }
  bool measuresLoc1() const {
    return (m_dirMask & (ProjectorMask::bendingDir)) != 0u;
  }
  bool measuresLoc0() const {
    return (m_dirMask & (ProjectorMask::nonBendingDir)) != 0u;
  }

  double driftRadius() const { return m_radius; }
  double time() const { return m_time; }
  const std::array<double, 3>& covariance() const { return m_cov; }

 private:
  Vector3 m_pos{Vector3::Zero()};
  Vector3 m_dir{Vector3::Zero()};
  Vector3 m_toNext{Vector3::Zero()};
  Vector3 m_planeNorm{m_dir.cross(m_toNext).normalized()};
  double m_radius{0.};
  double m_time{0.};
  ProjectorMask m_dirMask{ProjectorMask::bothDirections};
  std::array<double, 3> m_cov{Acts::filledArray<double, 3>(0.)};
  bool m_isStraw{true};
};

BOOST_AUTO_TEST_SUITE(StrawLineSeederTest)

void testResidual(const Pars_t& linePars, const TestSpacePoint& testPoint) {
  constexpr auto logLvl = Logging::Level::INFO;
  using namespace Acts::detail::LineHelper;
  using ResidualIdx = StrawLineFitAuxiliaries::ResidualIdx;
  Config_t resCfg{};
  resCfg.useHessian = true;
  resCfg.calcAlongStrip = false;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::y0, ParIdx::phi, ParIdx::theta};
  Line_t line{};
  line.updateParameters(linePars);

  std::cout << "\n\n\nResidual test - Test line: " << toString(line.position())
            << ", " << toString(line.direction()) << std::endl;

  StrawLineFitAuxiliaries resCalc{resCfg,
                                  Acts::getDefaultLogger("testRes", logLvl)};
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
              << ", to-next: " << toString(testPoint.toNextSensor())
              << std::endl;
    const Vector3& n{testPoint.planeNormal()};
    const Vector3 planeIsect = intersectPlane(line.position(), line.direction(),
                                              n, testPoint.localPosition())
                                   .position();
    const Vector3 delta = (planeIsect - testPoint.localPosition());

    std::cout << "Residual: " << toString(resCalc.residual())
              << ", delta: " << toString(delta)
              << ", <delta, n> =" << delta.dot(n) << std::endl;

    Acts::ActsSquareMatrix<3> coordTrf{Acts::ActsSquareMatrix<3>::Zero()};
    coordTrf.col(ResidualIdx::bending) = testPoint.measuresLoc1()
                                             ? testPoint.toNextSensor()
                                             : testPoint.sensorDirection();
    coordTrf.col(ResidualIdx::nonBending) = testPoint.measuresLoc1()
                                                ? testPoint.sensorDirection()
                                                : testPoint.toNextSensor();
    coordTrf.col(2) = n;
    std::cout << "Coordinate trf: \n" << coordTrf << std::endl;

    const Vector3 trfDelta = coordTrf.inverse() * delta;
    std::cout << "Transformed delta: " << toString(trfDelta) << std::endl;

    if (testPoint.measuresLoc1()) {
      BOOST_CHECK_CLOSE(trfDelta[ResidualIdx::bending],
                        resCalc.residual()[ResidualIdx::bending], 1.e-10);
    } else {
      BOOST_CHECK_CLOSE(trfDelta[ResidualIdx::bending] *
                            static_cast<double>(resCfg.calcAlongStrip),
                        resCalc.residual()[ResidualIdx::bending], 1.e-10);
    }

    if (testPoint.measuresLoc0()) {
      BOOST_CHECK_CLOSE(trfDelta[ResidualIdx::nonBending],
                        resCalc.residual()[ResidualIdx::nonBending], 1.e-10);
    } else {
      BOOST_CHECK_CLOSE(trfDelta[ResidualIdx::nonBending] *
                            static_cast<double>(resCfg.calcAlongStrip),
                        resCalc.residual()[ResidualIdx::nonBending], 1.e-10);
    }
  }

  constexpr double h = 5.e-9;
  constexpr double tolerance = 1.e-3;
  for (auto par : resCfg.parsToUse) {
    Pars_t lineParsUp{linePars}, lineParsDn{linePars};
    lineParsUp[static_cast<std::size_t>(par)] += h;
    lineParsDn[static_cast<std::size_t>(par)] -= h;
    Line_t lineUp{}, lineDn{};
    lineUp.updateParameters(lineParsUp);
    lineDn.updateParameters(lineParsDn);

    StrawLineFitAuxiliaries resCalcUp{
        resCfg, Acts::getDefaultLogger("testResUp", logLvl)};
    StrawLineFitAuxiliaries resCalcDn{
        resCfg, Acts::getDefaultLogger("testResDn", logLvl)};

    resCalcUp.updateSpatialResidual(lineUp, testPoint);
    resCalcDn.updateSpatialResidual(lineDn, testPoint);

    const Vector numDeriv =
        (resCalcUp.residual() - resCalcDn.residual()) / (2. * h);

    std::cout << "Derivative test: " << StrawLineFitAuxiliaries::parName(par)
              << ", derivative: " << toString(resCalc.gradient(par))
              << ",  numerical: " << toString(numDeriv) << std::endl;
    BOOST_CHECK_LE((numDeriv - resCalc.gradient(par)).norm() /
                       std::max(numDeriv.norm(), 1.),
                   tolerance);
    /// Next attempt to calculate the second derivative
    for (auto par1 : resCfg.parsToUse) {
      lineParsUp = linePars;
      lineParsDn = linePars;
      lineParsUp[static_cast<std::size_t>(par1)] += h;
      lineParsDn[static_cast<std::size_t>(par1)] -= h;

      lineUp.updateParameters(lineParsUp);
      lineDn.updateParameters(lineParsDn);

      resCalcUp.updateSpatialResidual(lineUp, testPoint);
      resCalcDn.updateSpatialResidual(lineDn, testPoint);

      const Vector numDeriv1{
          (resCalcUp.gradient(par) - resCalcDn.gradient(par)) / (2. * h)};
      const Vector& analyticDeriv = resCalc.hessian(par, par1);
      std::cout << "Second deriv (" << StrawLineFitAuxiliaries::parName(par)
                << ", " << StrawLineFitAuxiliaries::parName(par1)
                << ") -- numerical: " << toString(numDeriv1)
                << ", analytic: " << toString(analyticDeriv) << std::endl;
      BOOST_CHECK_LE((numDeriv1 - analyticDeriv).norm(), tolerance);
    }
  }
}

BOOST_AUTO_TEST_CASE(WireResidualTest) {
  /// Set the line to be 45 degrees
  using Pars_t = Line_t::ParamVector;
  using ParIdx = Line_t::ParIndex;
  Pars_t linePars{};
  linePars[static_cast<std::size_t>(ParIdx::phi)] = 90. * 1_degree;
  linePars[static_cast<std::size_t>(ParIdx::theta)] = 45. * 1_degree;

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

  linePars[static_cast<std::size_t>(ParIdx::phi)] = 30. * 1_degree;
  linePars[static_cast<std::size_t>(ParIdx::theta)] = 60. * 1_degree;

  testResidual(linePars,
               TestSpacePoint{Vector{100. * 1_cm, 50. * 1_cm, 30. * 1_cm},
                              wireDir, 10. * 1_cm});

  testResidual(linePars,
               TestSpacePoint{Vector{100. * 1_cm, 50. * 1_cm, 30. * 1_cm},
                              wireDir, 10. * 1_cm, true});

  linePars[static_cast<std::size_t>(ParIdx::phi)] = 60. * 1_degree;
  linePars[static_cast<std::size_t>(ParIdx::theta)] = 30. * 1_degree;

  testResidual(linePars,
               TestSpacePoint{Vector{100. * 1_cm, 50. * 1_cm, 30. * 1_cm},
                              wireDir, 10. * 1_cm});
  testResidual(linePars,
               TestSpacePoint{Vector{100. * 1_cm, 50. * 1_cm, 30. * 1_cm},
                              wireDir, 10. * 1_cm, true});

  testResidual(
      linePars,
      TestSpacePoint{Vector{100. * 1_cm, 50. * 1_cm, 30. * 1_cm},
                     makeDirectionFromPhiTheta(20. * 1_degree, 30 * 1_degree),
                     10. * 1_cm});
  testResidual(
      linePars,
      TestSpacePoint{Vector{100. * 1_cm, 50. * 1_cm, 30. * 1_cm},
                     makeDirectionFromPhiTheta(40. * 1_degree, 50 * 1_degree),
                     10. * 1_cm, true});
}

BOOST_AUTO_TEST_CASE(StripResidual) {
  Pars_t linePars{};
  linePars[static_cast<std::size_t>(ParIdx::phi)] = 60. * 1_degree;
  linePars[static_cast<std::size_t>(ParIdx::theta)] = 45 * 1_degree;

  testResidual(
      linePars,
      TestSpacePoint{Vector{75. * 1_cm, -75. * 1_cm, 100. * 1_cm},
                     makeDirectionFromPhiTheta(90. * 1_degree, 90. * 1_degree),
                     makeDirectionFromPhiTheta(0. * 1_degree, 90 * 1_degree),
                     TestSpacePoint::bendingDir});

  testResidual(
      linePars,
      TestSpacePoint{Vector{75. * 1_cm, -75. * 1_cm, 100. * 1_cm},
                     makeDirectionFromPhiTheta(0. * 1_degree, 90. * 1_degree),
                     makeDirectionFromPhiTheta(90. * 1_degree, 90 * 1_degree),
                     TestSpacePoint::nonBendingDir});

  testResidual(
      linePars,
      TestSpacePoint{Vector{75. * 1_cm, -75. * 1_cm, 100. * 1_cm},
                     makeDirectionFromPhiTheta(90. * 1_degree, 90. * 1_degree),
                     makeDirectionFromPhiTheta(0. * 1_degree, 90 * 1_degree),
                     TestSpacePoint::bothDirections});

  testResidual(
      linePars,
      TestSpacePoint{Vector{75. * 1_cm, -75. * 1_cm, 100. * 1_cm},
                     makeDirectionFromPhiTheta(30. * 1_degree, 90. * 1_degree),
                     makeDirectionFromPhiTheta(60. * 1_degree, 90 * 1_degree),
                     TestSpacePoint::bothDirections});
}

BOOST_AUTO_TEST_CASE(CombinatorialSeedSolverStripsTest) {
  RandomEngine rndEngine{23568};
  const std::size_t nStrips = 8;
  const std::size_t nEvents = 100;

  const std::array<Vector3, nStrips> stripDirections = {
      Vector3::UnitX(),
      Vector3::UnitX(),
      makeDirectionFromPhiTheta(1.5 * 1_degree, 90. * 1_degree),
      makeDirectionFromPhiTheta(-1.5 * 1_degree, 90. * 1_degree),
      makeDirectionFromPhiTheta(1.5 * 1_degree, 90. * 1_degree),
      makeDirectionFromPhiTheta(-1.5 * 1_degree, 90. * 1_degree),
      Vector3::UnitX(),
      Vector3::UnitX()};

  constexpr std::array<double, nStrips> distancesZ{-20., -264., 0.,  300.,
                                                   350,  370.,  400, 450.};

  std::array<double, nStrips> parameters{};
  std::array<Vector3, nStrips> intersections{};

  // pseudo track initilization
  Line_t line{};
  Pars_t linePars{};
  linePars[static_cast<std::size_t>(ParIdx::x0)] = 0. * 1_mm;
  linePars[static_cast<std::size_t>(ParIdx::y0)] = 0. * 1_mm;
  linePars[static_cast<std::size_t>(ParIdx::phi)] = 0. * 1_degree;
  linePars[static_cast<std::size_t>(ParIdx::theta)] = 0. * 1_degree;
  line.updateParameters(linePars);

  for (std::size_t i = 0; i < nEvents; i++) {
    std::cout << "Processing Event: " << i << std::endl;
    // update pseudo track parameters with random values
    linePars[static_cast<double>(ParIdx::x0)] = rndEngine() % 1000 - 500.;
    linePars[static_cast<double>(ParIdx::y0)] = rndEngine() % 1000 - 500.;
    linePars[static_cast<double>(ParIdx::phi)] = (rndEngine() % 90) * 1_degree;
    linePars[static_cast<double>(ParIdx::theta)] =
        (rndEngine() % 90) * 1_degree;
    line.updateParameters(linePars);

    Vector3 muonPos = line.position();
    Vector3 muonDir = line.direction();

    std::array<std::unique_ptr<TestSpacePoint>, nStrips> spacePoints{};
    for (std::size_t layer = 0; layer < nStrips; layer++) {
      auto intersection = PlanarHelper::intersectPlane(
          muonPos, muonDir, Vector3::UnitZ(), distancesZ[layer]);
      intersections[layer] = intersection.position();
      // where the hit is along the strip
      auto alongStrip = PlanarHelper::intersectPlane(
          intersections[layer], stripDirections[layer], Vector3::UnitX(), 0.);
      parameters[layer] = alongStrip.pathLength();

      Vector3 stripPos =
          intersections[layer] + parameters[layer] * stripDirections[layer];

      spacePoints[layer] = std::make_unique<TestSpacePoint>(
          stripPos, stripDirections[layer], Vector3::UnitZ(),
          TestSpacePoint::nonBendingDir);
    }

    // let's pick four strips on fours different layers
    // check combinatorics
    std::array<TestSpacePoint*, 4> seedSpacePoints = {};
    for (unsigned int l = 0; l < distancesZ.size() - 3; ++l) {
      seedSpacePoints[0] = spacePoints[l].get();
      for (unsigned int k = l + 1; k < distancesZ.size() - 2; ++k) {
        seedSpacePoints[1] = spacePoints[k].get();
        for (unsigned int m = k + 1; m < distancesZ.size() - 1; ++m) {
          seedSpacePoints[2] = spacePoints[m].get();
          for (unsigned int n = m + 1; n < distancesZ.size(); ++n) {
            seedSpacePoints[3] = spacePoints[n].get();

            const SquareMatrix2 bMatrix =
                CombinatorialSeedSolver::betaMatrix(seedSpacePoints);
            if (std::abs(bMatrix.determinant()) <
                std::numeric_limits<float>::epsilon()) {
              std::cout
                  << "Beta Matrix has zero determinant -skip the combination"
                  << std::endl;
              continue;
            }
            const std::array<double, 4> distsAlongStrip =
                CombinatorialSeedSolver::defineParameters(bMatrix,
                                                          seedSpacePoints);
            const auto [seedPos, seedDir] =
                CombinatorialSeedSolver::seedSolution(seedSpacePoints,
                                                      distsAlongStrip);

            // TO DO include the tests/ checks
          }
        }
      }
    }
  }
}

}  // namespace Acts::Test
BOOST_AUTO_TEST_SUITE_END()

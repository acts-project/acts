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
#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "Acts/Utilities/detail/Polynomials.hpp"

#include <random>

using namespace Acts;
using namespace Acts::detail;
using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;
using namespace Acts::UnitLiterals;
using namespace Acts::PlanarHelper;
using RandomEngine = std::mt19937;

using Line_t = CompSpacePointAuxiliaries::Line_t;
using Config_t = CompSpacePointAuxiliaries::Config;
using ParIdx = CompSpacePointAuxiliaries::FitParIndex;

using Vector = Line_t::Vector;
using Pars_t = Line_t::ParamVector;

constexpr auto logLvl = Logging::Level::VERBOSE;

namespace Acts::Test {
class TestSpacePoint {
 public:
  /// @brief Mask to express which track directions are measured
  ///        by the space point
  enum ProjectorMask : std::size_t {
    bendingDir = 1 << 0,
    nonBendingDir = 1 << 1,
    bothDirections = bendingDir | nonBendingDir
  };
  /// @brief Constructor for straw-like space points
  /// @param pos: A point on the straw wire
  /// @param wireDir: Orientation of the straw in space
  /// @param radius: Drift-radius of the straw measurement
  /// @param measNonPrec: Flag toggling whether the position along
  ///                     the wire shall be taken into account for
  ///                     residual calculation
  /// @param cov: Straw covariance
  TestSpacePoint(
      const Vector3& pos, const Vector3& wireDir, const double radius,
      bool measNonPrec = false,
      const std::array<double, 3>& cov = Acts::filledArray<double, 3>(0.))
      : m_pos{pos},
        m_dir{wireDir},
        m_radius{radius},
        m_dirMask{measNonPrec ? bothDirections : bendingDir},
        m_cov{cov} {}

  /// @brief Constructor for strip space points without time measurements
  /// @param pos: Position of the (combined) strip space point
  /// @param stripDir: Orientation of the strip in space
  /// @param toNext: Vector pointing to the next strip
  /// @param mask: Enum Toggling whether the strip measures loc0 or loc1
  /// @param cov: Strip's covariance
  TestSpacePoint(
      const Vector3& pos, const Vector3& stripDir, const Vector3& toNext,
      ProjectorMask mask,
      const std::array<double, 3>& cov = Acts::filledArray<double, 3>(0.))
      : m_pos{pos},
        m_dir{stripDir},
        m_toNext{toNext},
        m_dirMask{mask},
        m_cov{cov},
        m_isStraw{false} {}
  ///
  TestSpacePoint(
      const Vector3& pos, const Vector3& stripDir, const Vector3& toNext,
      const double stripTime,
      const std::array<double, 3>& cov = Acts::filledArray<double, 3>(0.))
      : m_pos{pos},
        m_dir{stripDir},
        m_toNext{toNext},
        m_time{stripTime},
        m_cov{cov},
        m_isStraw{false} {}
  const Vector3& localPosition() const { return m_pos; }
  const Vector3& sensorDirection() const { return m_dir; }
  const Vector3& toNextSensor() const { return m_toNext; }
  const Vector3& planeNormal() const { return m_planeNorm; }

  bool isStraw() const { return m_isStraw; }
  bool hasTime() const { return m_time.has_value(); }
  bool measuresLoc1() const {
    return (m_dirMask & (ProjectorMask::bendingDir)) != 0u;
  }
  bool measuresLoc0() const {
    return (m_dirMask & (ProjectorMask::nonBendingDir)) != 0u;
  }

  double driftRadius() const { return m_radius; }
  double time() const { return m_time.value_or(0.); }
  const std::array<double, 3>& covariance() const { return m_cov; }

 private:
  Vector3 m_pos{Vector3::Zero()};
  Vector3 m_dir{Vector3::Zero()};
  Vector3 m_toNext{Vector3::Zero()};
  Vector3 m_planeNorm{m_dir.cross(m_toNext).normalized()};
  double m_radius{0.};
  std::optional<double> m_time{std::nullopt};
  ProjectorMask m_dirMask{ProjectorMask::bothDirections};
  std::array<double, 3> m_cov{Acts::filledArray<double, 3>(0.)};
  bool m_isStraw{true};
};

BOOST_AUTO_TEST_SUITE(StrawLineSeederTest)

void testResidual(const Pars_t& linePars, const TestSpacePoint& testPoint) {
  using namespace Acts::detail::LineHelper;
  using ResidualIdx = CompSpacePointAuxiliaries::ResidualIdx;
  Config_t resCfg{};
  resCfg.useHessian = true;
  resCfg.calcAlongStrip = false;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::y0, ParIdx::phi, ParIdx::theta};
  Line_t line{};
  line.updateParameters(linePars);

  std::cout << "\n\n\nResidual test - Test line: " << toString(line.position())
            << ", " << toString(line.direction()) << std::endl;

  CompSpacePointAuxiliaries resCalc{resCfg,
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

    CompSpacePointAuxiliaries resCalcUp{
        resCfg, Acts::getDefaultLogger("testResUp", logLvl)};
    CompSpacePointAuxiliaries resCalcDn{
        resCfg, Acts::getDefaultLogger("testResDn", logLvl)};

    resCalcUp.updateSpatialResidual(lineUp, testPoint);
    resCalcDn.updateSpatialResidual(lineDn, testPoint);

    const Vector numDeriv =
        (resCalcUp.residual() - resCalcDn.residual()) / (2. * h);

    std::cout << "Derivative test: " << CompSpacePointAuxiliaries::parName(par)
              << ", derivative: " << toString(resCalc.gradient(par))
              << ",  numerical: " << toString(numDeriv) << std::endl;
    BOOST_CHECK_LE((numDeriv - resCalc.gradient(par)).norm() /
                       std::max(numDeriv.norm(), 1.),
                   tolerance);
    /// Next attempt to calculate the second derivative
    for (auto par1 : resCfg.parsToUse) {
      const Vector numDeriv1{
          (resCalcUp.gradient(par1) - resCalcDn.gradient(par1)) / (2. * h)};
      const Vector& analyticDeriv = resCalc.hessian(par, par1);
      std::cout << "Second deriv (" << CompSpacePointAuxiliaries::parName(par)
                << ", " << CompSpacePointAuxiliaries::parName(par1)
                << ") -- numerical: " << toString(numDeriv1)
                << ", analytic: " << toString(analyticDeriv) << std::endl;
      BOOST_CHECK_LE((numDeriv1 - analyticDeriv).norm(), tolerance);
    }
  }
}


void testSeed(const std::array<TestSpacePoint*, 4>& spacePoints,
              const std::array<double, 4>& truthDistances,
              const Vector3& truthPosZ0, const Vector3& truthDir) {
  const SquareMatrix2 bMatrix =
      CombinatorialSeedSolver::betaMatrix(spacePoints);
  if (std::abs(bMatrix.determinant()) < std::numeric_limits<float>::epsilon()) {
    std::cout << "Beta Matrix has zero determinant -skip the combination"
              << std::endl;
    return;
  }
  const std::array<double, 4> distsAlongStrip =
      CombinatorialSeedSolver::defineParameters(bMatrix, spacePoints);
  const auto [seedPos, seedDir] =
      CombinatorialSeedSolver::seedSolution(spacePoints, distsAlongStrip);

  std::cout << " Reconstructed position " << toString(seedPos)
            << " and direction " << toString(seedDir) << std::endl;
  std::cout << " Truth position " << toString(truthPosZ0)
            << " and truth direction " << toString(truthDir) << std::endl;

  // check the distances along the strips if they are the same
  for (std::size_t i = 0; i < truthDistances.size(); i++) {
    BOOST_CHECK_LE(distsAlongStrip[i] + truthDistances[i],
                   std ::numeric_limits<float>::epsilon());
  }

  BOOST_CHECK_LE((seedPos - truthPosZ0).norm(),
                 std::numeric_limits<float>::epsilon());

  // check the direction
  BOOST_CHECK_LE(std::abs(std::abs(seedDir.dot(truthDir)) - 1.),
                 std::numeric_limits<float>::epsilon());
}


void timeStripResidualTest(const Pars_t& linePars, const double timeT0,
                           const TestSpacePoint& sp,
                           const Acts::Transform3& locToGlob) {
  using namespace Acts::detail::LineHelper;
  using ResidualIdx = CompSpacePointAuxiliaries::ResidualIdx;
  Config_t resCfg{};
  resCfg.localToGlobal = locToGlob;
  resCfg.useHessian = true;
  resCfg.calcAlongStrip = true;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::y0, ParIdx::phi, ParIdx::theta,
                      ParIdx::t0};
  Line_t line{};
  line.updateParameters(linePars);

  std::cout << "\n\n\nResidual test for line: " << toString(line.position())
            << ", " << toString(line.direction()) << " with t0: " << timeT0
            << "." << std::endl;

  CompSpacePointAuxiliaries resCalc{resCfg,
                                    Acts::getDefaultLogger("timeRes", logLvl)};
  resCalc.updateFullResidual(line, timeT0, sp);

  const Vector3 planeIsect =
      intersectPlane(line.position(), line.direction(), sp.planeNormal(),
                     sp.localPosition())
          .position();
  const double ToF =
      (locToGlob * planeIsect).norm() / Acts::PhysicalConstants::c + timeT0;

  BOOST_CHECK_CLOSE(resCalc.residual()[ResidualIdx::time], sp.time() - ToF,
                    1.e-10);
  std::cout << "Time of flight: " << ToF << ", measured time : " << sp.time()
            << " --> residual: " << (sp.time() - ToF)
            << " vs. from calculator: " << resCalc.residual()[ResidualIdx::time]
            << "." << std::endl;
  constexpr double h = 5.e-9;
  constexpr double tolerance = 1.e-3;

  Line_t lineUp{}, lineDn{};
  for (const auto partial : resCfg.parsToUse) {
    CompSpacePointAuxiliaries resCalcUp{
        resCfg, Acts::getDefaultLogger("timeResUp", logLvl)};
    CompSpacePointAuxiliaries resCalcDn{
        resCfg, Acts::getDefaultLogger("timeResDn", logLvl)};

    Pars_t lineParsUp{linePars}, lineParsDn{linePars};

    if (partial != ParIdx::t0) {
      lineParsUp[static_cast<std::size_t>(partial)] += h;
      lineParsDn[static_cast<std::size_t>(partial)] -= h;
      lineUp.updateParameters(lineParsUp);
      lineDn.updateParameters(lineParsDn);
      resCalcUp.updateFullResidual(lineUp, timeT0, sp);
      resCalcDn.updateFullResidual(lineDn, timeT0, sp);
    } else {
      resCalcUp.updateFullResidual(line, timeT0 + h, sp);
      resCalcDn.updateFullResidual(line, timeT0 - h, sp);
    }

    const Vector numDeriv =
        (resCalcUp.residual() - resCalcDn.residual()) / (2. * h);

    std::cout << "Derivative test: "
              << CompSpacePointAuxiliaries::parName(partial)
              << ", derivative: " << toString(resCalc.gradient(partial))
              << ",  numerical: " << toString(numDeriv) << std::endl;
    BOOST_CHECK_LE((numDeriv - resCalc.gradient(partial)).norm() /
                       std::max(numDeriv.norm(), 1.),
                   tolerance);
    for (const auto partial2 : resCfg.parsToUse) {
      const Vector numDeriv1{
          (resCalcUp.gradient(partial2) - resCalcDn.gradient(partial2)) /
          (2. * h)};
      const Vector& analyticDeriv = resCalc.hessian(partial, partial2);
      std::cout << "Second deriv ("
                << CompSpacePointAuxiliaries::parName(partial) << ", "
                << CompSpacePointAuxiliaries::parName(partial2)
                << ") -- numerical: " << toString(numDeriv1)
                << ", analytic: " << toString(analyticDeriv) << std::endl;
      BOOST_CHECK_LE((numDeriv1 - analyticDeriv).norm(), tolerance);
    }
  }
}
BOOST_AUTO_TEST_CASE(StrawDriftTimeCase) {
  using namespace Acts::detail::LineHelper;

  std::cout << "\n\n\n\nCalibration straw point test " << std::endl;

  constexpr double recordedT = 15._ns;
  constexpr std::array<double, 3> rtCoeffs{0., 75. / 1_ns,
                                           5000. / Acts::pow(1_ns, 2)};

  Config_t resCfg{};
  resCfg.useHessian = true;
  resCfg.calcAlongStrip = true;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::phi, ParIdx::y0, ParIdx::theta,
                      ParIdx::t0};

  auto makeCalculator = [&rtCoeffs, &resCfg](
                            const Pars_t& linePars, const Vector& pos,
                            const Vector& dir, const double t0) {
    Line_t line{};
    line.updateParameters(linePars);
    std::cout << "Calculate residual w.r.t. " << toString(line.position())
              << ", " << toString(line.direction()) << std::endl;
    auto isectP = lineIntersect(pos, dir, line.position(), line.direction());

    const double driftT =
        recordedT - isectP.position().norm() / PhysicalConstants::c - t0;

    const double driftR = polynomialSum(driftT, rtCoeffs);
    const double driftV = derivativeSum<1>(driftT, rtCoeffs);
    const double driftA = derivativeSum<2>(driftT, rtCoeffs);
    std::cout << "Create new space point @ " << toString(pos) << ", "
              << toString(dir) << ", using t0: " << t0 / 1._ns
              << " --> drfitT: " << driftT / 1._ns << " --> driftR: " << driftR
              << ", driftV: " << driftV << ", "
              << "driftA: " << driftA << std::endl;

    CompSpacePointAuxiliaries resCalc{
        resCfg, Acts::getDefaultLogger("timeRes", logLvl)};

    resCalc.updateFullResidual(line, t0, TestSpacePoint{pos, dir, driftR},
                               driftV, driftA);
    return resCalc;
  };

  Pars_t linePars{};
  linePars[static_cast<std::size_t>(ParIdx::phi)] = 90._degree;
  linePars[static_cast<std::size_t>(ParIdx::theta)] = 45_degree;
  linePars[static_cast<std::size_t>(ParIdx::x0)] = 0._cm;
  linePars[static_cast<std::size_t>(ParIdx::y0)] = -105_cm;

  const Vector wPos{0._cm, -75._cm, 150._cm};
  const Vector wDir{Vector::UnitX()};
  const double t0{10._ns};
  auto resCalc = makeCalculator(linePars, wPos, wDir, t0);
  std::cout << "Residual: " << toString(resCalc.residual()) << std::endl;
  for (const auto partial : resCfg.parsToUse) {
    std::cout << " *** partial: " << CompSpacePointAuxiliaries::parName(partial)
              << " " << toString(resCalc.gradient(partial)) << std::endl;
  }
  constexpr double h = 1.e-7;
  constexpr double tolerance = 1.e-3;
  for (const auto partial : resCfg.parsToUse) {
    Pars_t lineParsUp{linePars}, lineParsDn{linePars};
    double t0Up{t0}, t0Dn{t0};
    if (partial != ParIdx::t0) {
      lineParsUp[static_cast<std::size_t>(partial)] += h;
      lineParsDn[static_cast<std::size_t>(partial)] -= h;
    } else {
      t0Up += h;
      t0Dn -= h;
    }
    auto resCalcUp = makeCalculator(lineParsUp, wPos, wDir, t0Up);
    auto resCalcDn = makeCalculator(lineParsDn, wPos, wDir, t0Dn);

    const Vector numDeriv =
        (resCalcUp.residual() - resCalcDn.residual()) / (2. * h);
    std::cout << "\nPartial " << CompSpacePointAuxiliaries::parName(partial)
              << " --  numerical: " << toString(numDeriv)
              << ", analytical: " << toString(resCalc.gradient(partial)) << "\n"
              << std::endl;
    BOOST_CHECK_LE((numDeriv - resCalc.gradient(partial)).norm() /
                       std::max(numDeriv.norm(), 1.),
                   tolerance);
    for (const auto partial2 : resCfg.parsToUse) {
      const Vector numDeriv1{
          (resCalcUp.gradient(partial2) - resCalcDn.gradient(partial2)) /
          (2. * h)};
      const Vector& analyticDeriv = resCalc.hessian(partial, partial2);
      std::cout << "Second deriv ("
                << CompSpacePointAuxiliaries::parName(partial) << ", "
                << CompSpacePointAuxiliaries::parName(partial2)
                << ") -- numerical: " << toString(numDeriv1)
                << ", analytic: " << toString(analyticDeriv) << std::endl;
      BOOST_CHECK_LE((numDeriv1 - analyticDeriv).norm(), tolerance);
    }
  }
}
BOOST_AUTO_TEST_CASE(WireResidualTest) {
  // Set the line to be 45 degrees
  using Pars_t = Line_t::ParamVector;
  using ParIdx = Line_t::ParIndex;
  Pars_t linePars{};
  linePars[static_cast<std::size_t>(ParIdx::phi)] = 30._degree;
  linePars[static_cast<std::size_t>(ParIdx::theta)] = 60._degree;

  // Generate the first test measurement
  testResidual(linePars, TestSpacePoint{Vector{100._cm, 50._cm, 30._cm},
                                        Vector::UnitX(), 10._cm});
  testResidual(linePars,
               TestSpacePoint{Vector{100._cm, 50._cm, 30._cm},
                              makeDirectionFromPhiTheta(10._degree, 30._degree),
                              10._cm});

  testResidual(linePars, TestSpacePoint{Vector{100._cm, 50._cm, 30._cm},
                                        Vector::UnitX(), 10._cm, true});

  linePars[static_cast<std::size_t>(ParIdx::phi)] = 30._degree;
  linePars[static_cast<std::size_t>(ParIdx::theta)] = 60._degree;

  testResidual(linePars, TestSpacePoint{Vector{100._cm, 50._cm, 30._cm},
                                        Vector::UnitX(), 10._cm});

  testResidual(linePars, TestSpacePoint{Vector{100._cm, 50._cm, 30._cm},
                                        Vector::UnitX(), 10._cm, true});

  linePars[static_cast<std::size_t>(ParIdx::phi)] = 60._degree;
  linePars[static_cast<std::size_t>(ParIdx::theta)] = 30._degree;

  testResidual(linePars, TestSpacePoint{Vector{100._cm, 50._cm, 30._cm},
                                        Vector::UnitX(), 10._cm});
  testResidual(linePars, TestSpacePoint{Vector{100._cm, 50._cm, 30._cm},
                                        Vector::UnitX(), 10._cm, true});

  testResidual(
      linePars,
      TestSpacePoint{Vector{100._cm, 50._cm, 30._cm},
                     makeDirectionFromPhiTheta(20._degree, 30_degree), 10._cm});
  testResidual(linePars,
               TestSpacePoint{Vector{100._cm, 50._cm, 30._cm},
                              makeDirectionFromPhiTheta(40._degree, 50_degree),
                              10._cm, true});
}

BOOST_AUTO_TEST_CASE(StripResidual) {
  Pars_t linePars{};
  linePars[static_cast<std::size_t>(ParIdx::phi)] = 60._degree;
  linePars[static_cast<std::size_t>(ParIdx::theta)] = 45_degree;

  testResidual(linePars,
               TestSpacePoint{Vector{75._cm, -75._cm, 100._cm},
                              makeDirectionFromPhiTheta(90._degree, 90._degree),
                              makeDirectionFromPhiTheta(0._degree, 90_degree),
                              TestSpacePoint::bendingDir});

  testResidual(linePars,
               TestSpacePoint{Vector{75._cm, -75._cm, 100._cm},
                              makeDirectionFromPhiTheta(0._degree, 90._degree),
                              makeDirectionFromPhiTheta(90._degree, 90_degree),
                              TestSpacePoint::nonBendingDir});

  testResidual(linePars,
               TestSpacePoint{Vector{75._cm, -75._cm, 100._cm},
                              makeDirectionFromPhiTheta(90._degree, 90._degree),
                              makeDirectionFromPhiTheta(0._degree, 90_degree),
                              TestSpacePoint::bothDirections});

  testResidual(linePars,
               TestSpacePoint{Vector{75._cm, -75._cm, 100._cm},
                              makeDirectionFromPhiTheta(30._degree, 90._degree),
                              makeDirectionFromPhiTheta(60._degree, 90_degree),
                              TestSpacePoint::bothDirections});
}

BOOST_AUTO_TEST_CASE(CombinatorialSeedSolverStripsTest) {
  RandomEngine rndEngine{23568};
  const std::size_t nStrips = 8;
  const std::size_t nEvents = 100;
  constexpr auto x0_idx = static_cast<std::size_t>(ParIdx::x0);
  constexpr auto y0_idx = static_cast<std::size_t>(ParIdx::y0);
  constexpr auto phi_idx = static_cast<std::size_t>(ParIdx::phi);
  constexpr auto theta_idx = static_cast<std::size_t>(ParIdx::theta);

  const std::array<Vector3, nStrips> stripDirections = {
      Vector3::UnitX(),
      Vector3::UnitX(),
      makeDirectionFromPhiTheta(1.5_degree, 90_degree),
      makeDirectionFromPhiTheta(-1.5_degree, 90_degree),
      makeDirectionFromPhiTheta(1.5_degree, 90_degree),
      makeDirectionFromPhiTheta(-1.5_degree, 90_degree),
      Vector3::UnitX(),
      Vector3::UnitX()};

  constexpr std::array<double, nStrips> distancesZ{-20., -264., 0.,  300.,
                                                   350,  370.,  400, 450.};

  std::array<double, nStrips> parameters{};
  std::array<Vector3, nStrips> intersections{};

  // pseudo track initialization
  Line_t line{};
  Pars_t linePars{};
  linePars[x0_idx] = 0. * 1_mm;
  linePars[y0_idx] = 0. * 1_mm;
  linePars[phi_idx] = 0. * 1_degree;
  linePars[theta_idx] = 0. * 1_degree;
  line.updateParameters(linePars);

  for (std::size_t i = 0; i < nEvents; i++) {
    std::cout << "\n\n\nCombinatorial Seed test - Processing Event: " << i
              << std::endl;
    // update pseudo track parameters with random values
    linePars[x0_idx] = rndEngine() % 1000 - 500.;
    linePars[y0_idx] = rndEngine() % 1000 - 500.;
    linePars[phi_idx] = (rndEngine() % 90) * 1_degree;
    linePars[theta_idx] = (rndEngine() % 90) * 1_degree;
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

            const std::array<double, 4> truthDistances = {
                parameters[l], parameters[k], parameters[m], parameters[n]};

            testSeed(seedSpacePoints, truthDistances, muonPos, muonDir);
          }
        }
      }
    }
  }
}

=======
BOOST_AUTO_TEST_CASE(TimeStripResidual) {
  Pars_t linePars{};
  linePars[static_cast<std::size_t>(ParIdx::phi)] = 60._degree;
  linePars[static_cast<std::size_t>(ParIdx::theta)] = 45_degree;

  Acts::Transform3 locToGlob{Acts::Transform3::Identity()};

  timeStripResidualTest(
      linePars, 10.,
      TestSpacePoint{Vector{75._cm, -75._cm, 100._cm},
                     makeDirectionFromPhiTheta(30._degree, 90._degree),
                     makeDirectionFromPhiTheta(60._degree, 90_degree), 15},
      locToGlob);

  timeStripResidualTest(
      linePars, 10.,
      TestSpacePoint{Vector{75._cm, -75_cm, 200_cm},
                     makeDirectionFromPhiTheta(0_degree, 90_degree),
                     makeDirectionFromPhiTheta(60._degree, 75_degree), 15},
      locToGlob);

  locToGlob.translation() = Vector{75._cm, -75._cm, -35._cm};
  timeStripResidualTest(
      linePars, 10.,
      TestSpacePoint{Vector{75._cm, -75._cm, 100._cm},
                     makeDirectionFromPhiTheta(30._degree, 90._degree),
                     makeDirectionFromPhiTheta(60._degree, 90_degree), 15},
      locToGlob);

  timeStripResidualTest(
      linePars, 10.,
      TestSpacePoint{Vector{75._cm, -75_cm, 200_cm},
                     makeDirectionFromPhiTheta(0_degree, 90_degree),
                     makeDirectionFromPhiTheta(60._degree, 75_degree), 15},
      locToGlob);
  locToGlob *= Acts::AngleAxis3{
      30_degree, makeDirectionFromPhiTheta(30_degree, -45_degree)};
  timeStripResidualTest(
      linePars, 10.,
      TestSpacePoint{Vector{75._cm, -75._cm, 100._cm},
                     makeDirectionFromPhiTheta(30._degree, 90._degree),
                     makeDirectionFromPhiTheta(60._degree, 90_degree), 15},
      locToGlob);

  timeStripResidualTest(
      linePars, 10.,
      TestSpacePoint{Vector{75._cm, -75_cm, 200_cm},
                     makeDirectionFromPhiTheta(0_degree, 90_degree),
                     makeDirectionFromPhiTheta(60._degree, 75_degree), 15},
      locToGlob);
}

}  // namespace Acts::Test
BOOST_AUTO_TEST_SUITE_END()

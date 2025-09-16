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

constexpr auto logLvl = Logging::Level::INFO;

ACTS_LOCAL_LOGGER(getDefaultLogger("StrawLineResidualTest", logLvl));

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
  /// @brief Constructor for strip space points with time measurements
  /// @param pos: Position of the (combined) strip space point
  /// @param stripDir: Orientation of the strip in space
  /// @param toNext: Vector pointing to the next strip
  /// @param stripTime: Recorded strip time used to calculate the residual
  /// @param cov: Strip's covariance
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

  friend std::ostream& operator<<(std::ostream& ostr,
                                  const TestSpacePoint& sp) {
    ostr << (sp.isStraw() ? "straw" : "strip") << " space point @ "
         << toString(sp.localPosition());
    if (sp.isStraw()) {
      ostr << " wire direction: " << toString(sp.sensorDirection())
           << ", drift Radius: " << sp.driftRadius();
    } else {
      ostr << ", sensor/toNext/normal: " << toString(sp.sensorDirection())
           << "/" << toString(sp.toNextSensor()) << "/"
           << toString(sp.planeNormal());
    }
    if (sp.hasTime()) {
      ostr << ", time: " << sp.time();
    }
    ostr << ", non-bending: " << (sp.measuresLoc0() ? "yay" : "nay");
    ostr << ", bending: " << (sp.measuresLoc1() ? "yay" : "nay");
    const auto& cov = sp.covariance();
    ostr << ", covariance: (" << cov[0] << ", " << cov[1] << ", " << cov[2]
         << "), ";
    return ostr;
  }

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
  Line_t line{linePars};

  ACTS_INFO(__func__ << "() - " << __LINE__
                     << ": Residual test w.r.t. test line: "
                     << toString(line.position()) << ", "
                     << toString(line.direction()));

  CompSpacePointAuxiliaries resCalc{resCfg,
                                    Acts::getDefaultLogger("testRes", logLvl)};
  resCalc.updateSpatialResidual(line, testPoint);
  ACTS_INFO(__func__ << "() - " << __LINE__ << ": Test residual w.r.t. "
                     << testPoint);
  if (testPoint.isStraw()) {
    const double lineDist =
        signedDistance(testPoint.localPosition(), testPoint.sensorDirection(),
                       line.position(), line.direction());

    ACTS_INFO(__func__ << "() - " << __LINE__
                       << ": Residual: " << toString(resCalc.residual())
                       << ", line distance: " << lineDist);
    ///
    BOOST_CHECK_CLOSE(resCalc.residual()[toUnderlying(ResidualIdx::bending)],
                      lineDist - testPoint.driftRadius(), 1.e-12);
  } else {
    const Vector3& n{testPoint.planeNormal()};
    const Vector3 planeIsect = intersectPlane(line.position(), line.direction(),
                                              n, testPoint.localPosition())
                                   .position();
    const Vector3 delta = (planeIsect - testPoint.localPosition());

    ACTS_INFO(__func__ << "() - " << __LINE__
                       << ": Residual: " << toString(resCalc.residual())
                       << ", delta: " << toString(delta)
                       << ", <delta, n> =" << delta.dot(n));

    Acts::ActsSquareMatrix<3> coordTrf{Acts::ActsSquareMatrix<3>::Zero()};
    coordTrf.col(toUnderlying(ResidualIdx::bending)) =
        testPoint.measuresLoc1() ? testPoint.toNextSensor()
                                 : testPoint.sensorDirection();
    coordTrf.col(toUnderlying(ResidualIdx::nonBending)) =
        testPoint.measuresLoc1() ? testPoint.sensorDirection()
                                 : testPoint.toNextSensor();
    coordTrf.col(2) = n;
    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Coordinate trf: \n"
                        << coordTrf);

    const Vector3 trfDelta = coordTrf.inverse() * delta;
    ACTS_DEBUG(__func__ << "() - " << __LINE__
                        << ": Transformed delta: " << toString(trfDelta));

    if (testPoint.measuresLoc1()) {
      BOOST_CHECK_CLOSE(trfDelta[toUnderlying(ResidualIdx::bending)],
                        resCalc.residual()[toUnderlying(ResidualIdx::bending)],
                        1.e-10);
    } else {
      BOOST_CHECK_CLOSE(trfDelta[toUnderlying(ResidualIdx::bending)] *
                            static_cast<double>(resCfg.calcAlongStrip),
                        resCalc.residual()[toUnderlying(ResidualIdx::bending)],
                        1.e-10);
    }
    if (testPoint.measuresLoc0()) {
      BOOST_CHECK_CLOSE(
          trfDelta[toUnderlying(ResidualIdx::nonBending)],
          resCalc.residual()[toUnderlying(ResidualIdx::nonBending)], 1.e-10);
    } else {
      BOOST_CHECK_CLOSE(
          trfDelta[toUnderlying(ResidualIdx::nonBending)] *
              static_cast<double>(resCfg.calcAlongStrip),
          resCalc.residual()[toUnderlying(ResidualIdx::nonBending)], 1.e-10);
    }
  }

  constexpr double h = 5.e-9;
  constexpr double tolerance = 1.e-3;
  for (auto par : resCfg.parsToUse) {
    Pars_t lineParsUp{linePars}, lineParsDn{linePars};
    lineParsUp[toUnderlying(par)] += h;
    lineParsDn[toUnderlying(par)] -= h;
    Line_t lineUp{lineParsUp}, lineDn{lineParsDn};

    CompSpacePointAuxiliaries resCalcUp{
        resCfg, Acts::getDefaultLogger("testResUp", logLvl)};
    CompSpacePointAuxiliaries resCalcDn{
        resCfg, Acts::getDefaultLogger("testResDn", logLvl)};

    resCalcUp.updateSpatialResidual(lineUp, testPoint);
    resCalcDn.updateSpatialResidual(lineDn, testPoint);

    const Vector numDeriv =
        (resCalcUp.residual() - resCalcDn.residual()) / (2. * h);

    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Derivative test: "
                        << CompSpacePointAuxiliaries::parName(par)
                        << ", derivative: " << toString(resCalc.gradient(par))
                        << ",  numerical: " << toString(numDeriv));
    BOOST_CHECK_LE((numDeriv - resCalc.gradient(par)).norm() /
                       std::max(numDeriv.norm(), 1.),
                   tolerance);
    /// Next attempt to calculate the second derivative
    for (auto par1 : resCfg.parsToUse) {
      const Vector numDeriv1{
          (resCalcUp.gradient(par1) - resCalcDn.gradient(par1)) / (2. * h)};
      const Vector& analyticDeriv = resCalc.hessian(par, par1);
      ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Second deriv ("
                          << CompSpacePointAuxiliaries::parName(par) << ", "
                          << CompSpacePointAuxiliaries::parName(par1)
                          << ") -- numerical: " << toString(numDeriv1)
                          << ", analytic: " << toString(analyticDeriv));
      BOOST_CHECK_LE((numDeriv1 - analyticDeriv).norm(), tolerance);
    }
  }
}

void testSeed(const std::array<std::shared_ptr<TestSpacePoint>, 4>& spacePoints,
              const std::array<double, 4>& truthDistances,
              const Vector3& truthPosZ0, const Vector3& truthDir) {
  const SquareMatrix2 bMatrix =
      CombinatorialSeedSolver::betaMatrix(spacePoints);
  if (std::abs(bMatrix.determinant()) < std::numeric_limits<float>::epsilon()) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": Beta Matrix \n"
                          << bMatrix
                          << "\nhas zero determinant - skip the combination");
    return;
  }
  const std::array<double, 4> distsAlongStrip =
      CombinatorialSeedSolver::defineParameters(bMatrix, spacePoints);
  const auto [seedPos, seedDir] =
      CombinatorialSeedSolver::seedSolution(spacePoints, distsAlongStrip);

  ACTS_DEBUG(__func__ << "() " << __LINE__ << ": Reconstructed position "
                      << toString(seedPos) << " and direction "
                      << toString(seedDir));
  ACTS_DEBUG(__func__ << "() " << __LINE__ << ": Truth position "
                      << toString(truthPosZ0) << " and truth direction "
                      << toString(truthDir));

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
  Line_t line{linePars};

  ACTS_INFO(__func__ << "() - " << __LINE__ << ": Residual test for line: "
                     << toString(line.position()) << ", "
                     << toString(line.direction()) << " with t0: " << timeT0
                     << ".");

  CompSpacePointAuxiliaries resCalc{resCfg,
                                    Acts::getDefaultLogger("timeRes", logLvl)};
  resCalc.updateFullResidual(line, timeT0, sp);

  const Vector3 planeIsect =
      intersectPlane(line.position(), line.direction(), sp.planeNormal(),
                     sp.localPosition())
          .position();
  const double ToF =
      (locToGlob * planeIsect).norm() / Acts::PhysicalConstants::c + timeT0;

  BOOST_CHECK_CLOSE(resCalc.residual()[toUnderlying(ResidualIdx::time)],
                    sp.time() - ToF, 1.e-10);
  ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Time of flight: " << ToF
                      << ", measured time : " << sp.time() << " --> residual: "
                      << (sp.time() - ToF) << " vs. from calculator: "
                      << resCalc.residual()[toUnderlying(ResidualIdx::time)]
                      << ".");
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
      lineParsUp[toUnderlying(partial)] += h;
      lineParsDn[toUnderlying(partial)] -= h;
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

    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Derivative test: "
                        << CompSpacePointAuxiliaries::parName(partial)
                        << ", derivative: "
                        << toString(resCalc.gradient(partial))
                        << ",  numerical: " << toString(numDeriv));
    BOOST_CHECK_LE((numDeriv - resCalc.gradient(partial)).norm() /
                       std::max(numDeriv.norm(), 1.),
                   tolerance);
    for (const auto partial2 : resCfg.parsToUse) {
      const Vector numDeriv1{
          (resCalcUp.gradient(partial2) - resCalcDn.gradient(partial2)) /
          (2. * h)};
      const Vector& analyticDeriv = resCalc.hessian(partial, partial2);
      ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Second deriv ("
                          << CompSpacePointAuxiliaries::parName(partial) << ", "
                          << CompSpacePointAuxiliaries::parName(partial2)
                          << ") -- numerical: " << toString(numDeriv1)
                          << ", analytic: " << toString(analyticDeriv));
      BOOST_CHECK_LE((numDeriv1 - analyticDeriv).norm(), tolerance);
    }
  }
}
BOOST_AUTO_TEST_CASE(StrawDriftTimeCase) {
  using namespace Acts::detail::LineHelper;

  ACTS_INFO("Calibration straw point test ");

  constexpr double recordedT = 15._ns;
  constexpr std::array<double, 3> rtCoeffs{0., 75. / 1_ns,
                                           5000. / Acts::pow(1_ns, 2)};

  Config_t resCfg{};
  resCfg.useHessian = true;
  resCfg.calcAlongStrip = true;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::phi, ParIdx::y0, ParIdx::theta,
                      ParIdx::t0};

  auto makeCalculator = [&rtCoeffs, &resCfg](
                            const std::string& calcName, const Pars_t& linePars,
                            const Vector& pos, const Vector& dir,
                            const double t0) {
    Line_t line{linePars};
    ACTS_INFO(__func__ << "() " << __LINE__ << ": Calculate residual w.r.t. "
                       << toString(line.position()) << ", "
                       << toString(line.direction()));
    auto isectP = lineIntersect(pos, dir, line.position(), line.direction());

    const double driftT = recordedT -
                          (resCfg.localToGlobal * isectP.position()).norm() /
                              PhysicalConstants::c -
                          t0;

    const double driftR = polynomialSum(driftT, rtCoeffs);
    const double driftV = derivativeSum<1>(driftT, rtCoeffs);
    const double driftA = derivativeSum<2>(driftT, rtCoeffs);
    ACTS_INFO(__func__ << "() " << __LINE__ << ": Create new space point @ "
                       << toString(pos) << ", " << toString(dir)
                       << ", using t0: " << t0 / 1._ns << " --> drfitT: "
                       << driftT / 1._ns << " --> driftR: " << driftR
                       << ", driftV: " << driftV << ", "
                       << "driftA: " << driftA);

    CompSpacePointAuxiliaries resCalc{resCfg,
                                      Acts::getDefaultLogger(calcName, logLvl)};

    resCalc.updateFullResidual(line, t0, TestSpacePoint{pos, dir, driftR},
                               driftV, driftA);
    return resCalc;
  };

  auto testTimingResidual = [&makeCalculator, &resCfg](
                                const Pars_t& linePars, const Vector& wPos,
                                const Vector& wDir, const double t0) {
    auto resCalc = makeCalculator("strawT0Res", linePars, wPos, wDir, t0);
    ACTS_DEBUG(__func__ << "() - " << __LINE__
                        << ": Residual: " << toString(resCalc.residual()));
    for (const auto partial : resCfg.parsToUse) {
      ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": *** partial: "
                          << CompSpacePointAuxiliaries::parName(partial) << " "
                          << toString(resCalc.gradient(partial)));
    }
    constexpr double h = 1.e-7;
    constexpr double tolerance = 1.e-3;
    for (const auto partial : resCfg.parsToUse) {
      Pars_t lineParsUp{linePars}, lineParsDn{linePars};
      double t0Up{t0}, t0Dn{t0};
      if (partial != ParIdx::t0) {
        lineParsUp[toUnderlying(partial)] += h;
        lineParsDn[toUnderlying(partial)] -= h;
      } else {
        t0Up += h;
        t0Dn -= h;
      }
      auto resCalcUp =
          makeCalculator("strawT0ResUp", lineParsUp, wPos, wDir, t0Up);
      auto resCalcDn =
          makeCalculator("strawT0ResDn", lineParsDn, wPos, wDir, t0Dn);

      const Vector numDeriv =
          (resCalcUp.residual() - resCalcDn.residual()) / (2. * h);
      ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Partial "
                          << CompSpacePointAuxiliaries::parName(partial)
                          << " --  numerical: " << toString(numDeriv)
                          << ", analytical: "
                          << toString(resCalc.gradient(partial)));

      BOOST_CHECK_LE((numDeriv - resCalc.gradient(partial)).norm() /
                         std::max(numDeriv.norm(), 1.),
                     tolerance);
      for (const auto partial2 : resCfg.parsToUse) {
        const Vector numDeriv1{
            (resCalcUp.gradient(partial2) - resCalcDn.gradient(partial2)) /
            (2. * h)};
        const Vector& analyticDeriv = resCalc.hessian(partial, partial2);
        ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Second deriv ("
                            << CompSpacePointAuxiliaries::parName(partial)
                            << ", "
                            << CompSpacePointAuxiliaries::parName(partial2)
                            << ") -- numerical: " << toString(numDeriv1)
                            << ", analytic: " << toString(analyticDeriv));
        BOOST_CHECK_LE((numDeriv1 - analyticDeriv).norm(), tolerance);
      }
    }
  };
  Pars_t linePars{};
  linePars[toUnderlying(ParIdx::phi)] = 90._degree;
  linePars[toUnderlying(ParIdx::theta)] = 45_degree;
  linePars[toUnderlying(ParIdx::x0)] = 0._cm;
  linePars[toUnderlying(ParIdx::y0)] = -105_cm;

  testTimingResidual(linePars, Vector{0._cm, -75._cm, 150._cm}, Vector::UnitX(),
                     10._ns);
  testTimingResidual(linePars, Vector{0._cm, -75._cm, 150._cm},
                     Vector{1., 1., 0.}.normalized(), 10._ns);
  linePars[toUnderlying(ParIdx::phi)] = 60._degree;
  testTimingResidual(linePars, Vector{0._cm, -75._cm, 150._cm}, Vector::UnitX(),
                     10._ns);
  testTimingResidual(linePars, Vector{0._cm, -75._cm, 150._cm},
                     Vector{1., 1., 0.}.normalized(), 10._ns);

  resCfg.localToGlobal.translation() = Vector{10._cm, 20._cm, -50._cm};
  testTimingResidual(linePars, Vector{0._cm, -75._cm, 150._cm}, Vector::UnitX(),
                     10._ns);
  testTimingResidual(linePars, Vector{0._cm, -75._cm, 150._cm},
                     Vector{1., 1., 0.}.normalized(), 10._ns);

  //// Next test the displacement
  resCfg.localToGlobal *= Acts::AngleAxis3{
      30_degree, makeDirectionFromPhiTheta(30_degree, -45_degree)};
  testTimingResidual(linePars, Vector{0._cm, -75._cm, 150._cm}, Vector::UnitX(),
                     10._ns);
  testTimingResidual(linePars, Vector{0._cm, -75._cm, 150._cm},
                     Vector{1., 1., 0.}.normalized(), 10._ns);
}
BOOST_AUTO_TEST_CASE(WireResidualTest) {
  // Set the line to be 45 degrees
  using Pars_t = Line_t::ParamVector;
  using ParIdx = Line_t::ParIndex;
  Pars_t linePars{};
  linePars[toUnderlying(ParIdx::phi)] = 30._degree;
  linePars[toUnderlying(ParIdx::theta)] = 60._degree;

  // Generate the first test measurement
  testResidual(linePars, TestSpacePoint{Vector{100._cm, 50._cm, 30._cm},
                                        Vector::UnitX(), 10._cm});
  testResidual(linePars,
               TestSpacePoint{Vector{100._cm, 50._cm, 30._cm},
                              makeDirectionFromPhiTheta(10._degree, 30._degree),
                              10._cm});

  testResidual(linePars, TestSpacePoint{Vector{100._cm, 50._cm, 30._cm},
                                        Vector::UnitX(), 10._cm, true});

  linePars[toUnderlying(ParIdx::phi)] = 30._degree;
  linePars[toUnderlying(ParIdx::theta)] = 60._degree;

  testResidual(linePars, TestSpacePoint{Vector{100._cm, 50._cm, 30._cm},
                                        Vector::UnitX(), 10._cm});

  testResidual(linePars, TestSpacePoint{Vector{100._cm, 50._cm, 30._cm},
                                        Vector::UnitX(), 10._cm, true});

  linePars[toUnderlying(ParIdx::phi)] = 60._degree;
  linePars[toUnderlying(ParIdx::theta)] = 30._degree;

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
  linePars[toUnderlying(ParIdx::phi)] = 60._degree;
  linePars[toUnderlying(ParIdx::theta)] = 45_degree;

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

BOOST_AUTO_TEST_CASE(TimeStripResidual) {
  Pars_t linePars{};
  linePars[toUnderlying(ParIdx::phi)] = 60._degree;
  linePars[toUnderlying(ParIdx::theta)] = 45_degree;

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

BOOST_AUTO_TEST_CASE(ChiSqEvaluation) {
  Pars_t linePars{};
  linePars[toUnderlying(ParIdx::phi)] = 30._degree;
  linePars[toUnderlying(ParIdx::theta)] = 75_degree;
  linePars[toUnderlying(ParIdx::x0)] = 10._cm;
  linePars[toUnderlying(ParIdx::y0)] = -10_cm;

  Config_t resCfg{};
  resCfg.useHessian = true;
  resCfg.calcAlongStrip = true;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::phi, ParIdx::y0, ParIdx::theta,
                      ParIdx::t0};

  Line_t line{linePars};

  CompSpacePointAuxiliaries resCalc{resCfg,
                                    Acts::getDefaultLogger("testRes", logLvl)};

  const double t0{1._ns};

  auto testChi2 = [&line, &t0, &resCalc,
                   &resCfg](const TestSpacePoint& testMe) {
    using ChiSq_t = CompSpacePointAuxiliaries::ChiSqWithDerivatives;
    ChiSq_t chi2{};
    ACTS_INFO(__func__ << "() " << __LINE__ << ": Test space point " << testMe);

    resCalc.updateFullResidual(line, t0, testMe);

    resCalc.updateChiSq(chi2, testMe.covariance());
    ACTS_DEBUG(__func__ << "() - " << __LINE__
                        << ": Calculated chi2: " << chi2.chi2 << ", Gradient: "
                        << toString(chi2.gradient) << ", Hessian: \n"
                        << chi2.hessian);
    BOOST_CHECK_CLOSE(CompSpacePointAuxiliaries::chi2Term(
                          line, resCfg.localToGlobal, t0, testMe),
                      chi2.chi2, 1.e-7);

    constexpr double h = 1.e-7;
    constexpr double tolerance = 1.e-3;
    for (const auto par : resCfg.parsToUse) {
      Pars_t lineParsUp{line.parameters()};
      Pars_t lineParsDn{line.parameters()};

      const auto dIdx = toUnderlying(par);
      double t0Up{t0}, t0Dn{t0};
      if (par != ParIdx::t0) {
        lineParsUp[dIdx] += h;
        lineParsDn[dIdx] -= h;
      } else {
        t0Up += h;
        t0Dn -= h;
      }
      ChiSq_t chi2Up{}, chi2Dn{};
      /// Calculate the updated chi2
      line.updateParameters(lineParsUp);
      resCalc.updateFullResidual(line, t0Up, testMe);
      resCalc.updateChiSq(chi2Up, testMe.covariance());

      line.updateParameters(lineParsDn);
      resCalc.updateFullResidual(line, t0Dn, testMe);
      resCalc.updateChiSq(chi2Dn, testMe.covariance());

      const double anaDeriv = chi2.gradient[dIdx];
      const double numDeriv = (chi2Up.chi2 - chi2Dn.chi2) / (2. * h);

      ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Partial "
                          << CompSpacePointAuxiliaries::parName(par)
                          << " --  numerical: " << numDeriv
                          << ", analytical: " << anaDeriv);
      BOOST_CHECK_LE(std::abs(numDeriv - anaDeriv) / std::max(numDeriv, 1.),
                     tolerance);
      for (const auto par2 : resCfg.parsToUse) {
        if (par2 > par) {
          break;
        }
        const auto dIdx2 = toUnderlying(par2);
        const double anaHess = chi2.hessian(dIdx, dIdx2);
        const double numHess =
            (chi2Up.gradient[dIdx2] - chi2Dn.gradient[dIdx2]) / (2. * h);

        ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Second deriv ("
                            << CompSpacePointAuxiliaries::parName(par) << ", "
                            << CompSpacePointAuxiliaries::parName(par2)
                            << ") -- numerical: " << numHess
                            << ", analytic: " << anaHess);
        BOOST_CHECK_LE(std::abs(anaHess - numHess), tolerance);
      }
    }
  };

  /// Test orthogonal strips
  testChi2(TestSpacePoint{
      line.point(20._cm) + 5._cm * line.direction().cross(Vector::UnitX()),
      makeDirectionFromPhiTheta(0_degree, 90._degree),
      makeDirectionFromPhiTheta(90._degree, 0._degree),
      15._ns,
      {Acts::pow(5._cm, 2), Acts::pow(10._cm, 2), Acts::pow(1._ns, 2)}});

  /// Test strips with stereo
  testChi2(TestSpacePoint{
      line.point(20._cm) + 5._cm * line.direction().cross(Vector::UnitX()),
      makeDirectionFromPhiTheta(0_degree, 45._degree),
      makeDirectionFromPhiTheta(60._degree, 0._degree),
      15._ns,
      {Acts::pow(5._cm, 2), Acts::pow(10._cm, 2), Acts::pow(1._ns, 2)}});
  //// Test ordinary straws
  testChi2(TestSpacePoint{
      line.point(20._cm) + 5._cm * line.direction().cross(Vector::UnitX()),
      Vector3::UnitX(),
      5._cm,
      false,
      {0., Acts::pow(10._mm, 2), 0.}});

  /// Test straws with information on the position along the tube
  testChi2(TestSpacePoint{line.point(20._cm) +
                              5._cm * line.direction().cross(Vector::UnitX()) +
                              4._cm * Vector::UnitX(),
                          Vector3::UnitX(),
                          5._cm,
                          true,
                          {Acts::pow(0.5_cm, 2), Acts::pow(10._mm, 2), 0.}});
}

BOOST_AUTO_TEST_CASE(CombinatorialSeedSolverStripsTest) {
  RandomEngine rndEngine{23568};
  const std::size_t nStrips = 8;
  const std::size_t nEvents = 100;
  constexpr auto x0_idx = toUnderlying(ParIdx::x0);
  constexpr auto y0_idx = toUnderlying(ParIdx::y0);
  constexpr auto phi_idx = toUnderlying(ParIdx::phi);
  constexpr auto theta_idx = toUnderlying(ParIdx::theta);

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

  for (std::size_t i = 0; i < nEvents; ++i) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__
                          << " - Combinatorial Seed test - Processing Event: "
                          << i);
    // update pseudo track parameters with random values
    Pars_t linePars{};
    linePars[x0_idx] = std::uniform_real_distribution{-500., 500.}(rndEngine);
    linePars[y0_idx] = std::uniform_real_distribution{-500., 500.}(rndEngine);
    linePars[phi_idx] =
        std::uniform_real_distribution{0.1_degree, 179.9_degree}(rndEngine);
    linePars[theta_idx] =
        std::uniform_real_distribution{-180_degree, 180_degree}(rndEngine);
    Line_t line{linePars};

    Vector3 muonPos = line.position();
    Vector3 muonDir = line.direction();

    std::array<std::shared_ptr<TestSpacePoint>, nStrips> spacePoints{};
    for (std::size_t layer = 0; layer < nStrips; ++layer) {
      auto intersection = PlanarHelper::intersectPlane(
          muonPos, muonDir, Vector3::UnitZ(), distancesZ[layer]);
      intersections[layer] = intersection.position();
      // where the hit is along the strip
      auto alongStrip = PlanarHelper::intersectPlane(
          intersections[layer], stripDirections[layer], Vector3::UnitX(), 0.);

      parameters[layer] = alongStrip.pathLength();

      Vector3 stripPos =
          intersections[layer] + parameters[layer] * stripDirections[layer];

      spacePoints[layer] = std::make_shared<TestSpacePoint>(
          stripPos, stripDirections[layer], Vector3::UnitZ(),
          TestSpacePoint::nonBendingDir);
    }

    // let's pick four strips on fours different layers
    // check combinatorics
    std::array<std::shared_ptr<TestSpacePoint>, 4> seedSpacePoints{};
    for (unsigned int l = 0; l < distancesZ.size() - 3; ++l) {
      seedSpacePoints[0] = spacePoints[l];
      for (unsigned int k = l + 1; k < distancesZ.size() - 2; ++k) {
        seedSpacePoints[1] = spacePoints[k];
        for (unsigned int m = k + 1; m < distancesZ.size() - 1; ++m) {
          seedSpacePoints[2] = spacePoints[m];
          for (unsigned int n = m + 1; n < distancesZ.size(); ++n) {
            seedSpacePoints[3] = spacePoints[n];

            const std::array<double, 4> truthDistances = {
                parameters[l], parameters[k], parameters[m], parameters[n]};

            testSeed(seedSpacePoints, truthDistances, muonPos, muonDir);
          }
        }
      }
    }
  }
}

}  // namespace Acts::Test
BOOST_AUTO_TEST_SUITE_END()

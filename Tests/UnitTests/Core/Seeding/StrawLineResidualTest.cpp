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
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/detail/Polynomials.hpp"

#include <random>

using namespace Acts;
using namespace Acts::detail;
using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;
using namespace Acts::UnitLiterals;
using namespace Acts::PlanarHelper;
using namespace Acts::detail::LineHelper;
using namespace Acts::VectorHelpers;
using RandomEngine = std::mt19937;
using uniform = std::uniform_real_distribution<double>;

using Line_t = CompSpacePointAuxiliaries::Line_t;
using Config_t = CompSpacePointAuxiliaries::Config;
using ParIdx = CompSpacePointAuxiliaries::FitParIndex;
using ResidualIdx = CompSpacePointAuxiliaries::ResidualIdx;

using Vector = Line_t::Vector;
using Pars_t = Line_t::ParamVector;

constexpr auto logLvl = Logging::Level::INFO;
constexpr std::size_t nEvents = 100;
constexpr double tolerance = 1.e-3;

ACTS_LOCAL_LOGGER(getDefaultLogger("StrawLineResidualTest", logLvl));

namespace ActsTests {

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

/// @brief Generate a random straight track with a flat distribution in theta & y0
///        Range in theta is [0, 180] degrees in steps of 0.1 excluding the
///        vicinity around 90 degrees.
///          In y0, the generated range is [-500, 500] mm
Line_t::ParamVector generateLine(RandomEngine& engine) {
  using ParIndex = Line_t::ParIndex;
  Line_t::ParamVector linePars{};
  linePars[toUnderlying(ParIndex::phi)] =
      uniform{-180._degree, 180._degree}(engine);
  linePars[toUnderlying(ParIndex::y0)] = uniform{-5000., 5000.}(engine);
  linePars[toUnderlying(ParIndex::x0)] = uniform{-5000., 5000.}(engine);
  linePars[toUnderlying(ParIndex::theta)] =
      uniform{0.1_degree, 179.9_degree}(engine);
  if (Acts::abs(linePars[toUnderlying(ParIndex::theta)] - 90._degree) <
      0.2_degree) {
    return generateLine(engine);
  }
  ACTS_DEBUG("Generated parameters -- "
             << "theta: "
             << (linePars[toUnderlying(ParIndex::theta)] / 1._degree)
             << ", phi: " << (linePars[toUnderlying(ParIndex::phi)] / 1._degree)
             << ", y0: " << linePars[toUnderlying(ParIndex::y0)]
             << ", x0: " << linePars[toUnderlying(ParIndex::x0)]);
  return linePars;
}

template <typename T>
constexpr T absMax(const T& a, const T& b) {
  return std::max(Acts::abs(a), Acts::abs(b));
}
template <typename T, typename... argsT>
constexpr T absMax(const T& a, const argsT... args) {
  return std::max(Acts::abs(a), absMax(args...));
}
double angle(const Vector& v1, const Vector& v2) {
  const double dots = Acts::abs(v1.dot(v2));
  return std::acos(std::clamp(dots, 0., 1.));
}

#define CHECK_CLOSE(a, b, tol) \
  BOOST_CHECK_LE(std::abs(a - b) / absMax(a, b, 1.), tol);

#define COMPARE_VECTORS(v1, v2, tol)                             \
  {                                                              \
    CHECK_CLOSE(v1[toUnderlying(ResidualIdx::bending)],          \
                v2[toUnderlying(ResidualIdx::bending)], tol);    \
    CHECK_CLOSE(v1[toUnderlying(ResidualIdx::nonBending)],       \
                v2[toUnderlying(ResidualIdx::nonBending)], tol); \
    CHECK_CLOSE(v1[toUnderlying(ResidualIdx::time)],             \
                v2[toUnderlying(ResidualIdx::time)], tol);       \
  }

BOOST_AUTO_TEST_SUITE(SeedingSuite)

void testResidual(const Pars_t& linePars, const TestSpacePoint& testPoint) {
  Config_t resCfg{};
  resCfg.useHessian = true;
  resCfg.calcAlongStrip = false;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::y0, ParIdx::phi, ParIdx::theta};
  Line_t line{linePars};

  ACTS_INFO(__func__ << "() - " << __LINE__
                     << ":\n##############################################\n###"
                        "###########################################\n"
                     << "Residual test w.r.t. test line: "
                     << toString(line.position()) << ", "
                     << toString(line.direction())
                     << ", theta: " << (theta(line.direction()) / 1._degree)
                     << ", phi: " << (phi(line.direction())) / 1._degree);

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
                       << ", line distance: " << lineDist
                       << ",  analog residual: "
                       << (lineDist - testPoint.driftRadius()));
    ///
    CHECK_CLOSE(resCalc.residual()[toUnderlying(ResidualIdx::bending)],
                (lineDist - testPoint.driftRadius()), 1.e-8);
    if (testPoint.measuresLoc0()) {
      auto closePoint =
          lineIntersect(line.position(), line.direction(),
                        testPoint.localPosition(), testPoint.sensorDirection());

      ACTS_INFO(__func__
                << "() - " << __LINE__ << ": Distance along wire "
                << closePoint.pathLength() << ", residual: "
                << resCalc.residual()[toUnderlying(ResidualIdx::nonBending)]);

      CHECK_CLOSE(resCalc.residual()[toUnderlying(ResidualIdx::nonBending)],
                  closePoint.pathLength(), 1.e-10);
    }
  } else {
    const Vector3& n{testPoint.planeNormal()};
    const Vector3 planeIsect = intersectPlane(line.position(), line.direction(),
                                              n, testPoint.localPosition())
                                   .position();
    const Vector3 delta = (planeIsect - testPoint.localPosition());

    ACTS_DEBUG(__func__ << "() - " << __LINE__
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
      CHECK_CLOSE(trfDelta[toUnderlying(ResidualIdx::bending)],
                  resCalc.residual()[toUnderlying(ResidualIdx::bending)],
                  1.e-10);
    } else {
      CHECK_CLOSE(trfDelta[toUnderlying(ResidualIdx::bending)] *
                      static_cast<double>(resCfg.calcAlongStrip),
                  resCalc.residual()[toUnderlying(ResidualIdx::bending)],
                  1.e-10);
    }
    if (testPoint.measuresLoc0()) {
      CHECK_CLOSE(trfDelta[toUnderlying(ResidualIdx::nonBending)],
                  resCalc.residual()[toUnderlying(ResidualIdx::nonBending)],
                  1.e-10);
    } else {
      CHECK_CLOSE(trfDelta[toUnderlying(ResidualIdx::nonBending)] *
                      static_cast<double>(resCfg.calcAlongStrip),
                  resCalc.residual()[toUnderlying(ResidualIdx::nonBending)],
                  1.e-10);
    }
  }

  constexpr double h = 2.e-7;
  for (auto par : resCfg.parsToUse) {
    Pars_t lineParsUp{linePars};
    Pars_t lineParsDn{linePars};
    lineParsUp[toUnderlying(par)] += h;
    lineParsDn[toUnderlying(par)] -= h;
    const Line_t lineUp{lineParsUp};
    const Line_t lineDn{lineParsDn};

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

    COMPARE_VECTORS(numDeriv, resCalc.gradient(par), tolerance);
    /// Next attempt to calculate the second derivative
    for (auto par1 : resCfg.parsToUse) {
      if (par1 > par) {
        break;
      }
      const Vector numDeriv1{
          (resCalcUp.gradient(par1) - resCalcDn.gradient(par1)) / (2. * h)};
      ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Second deriv ("
                          << CompSpacePointAuxiliaries::parName(par) << ", "
                          << CompSpacePointAuxiliaries::parName(par1)
                          << ") -- numerical: " << toString(numDeriv1)
                          << ", analytic: "
                          << toString(resCalc.hessian(par, par1)));
      COMPARE_VECTORS(numDeriv1, resCalc.hessian(par, par1), tolerance);
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
  for (std::size_t i = 0; i < truthDistances.size(); ++i) {
    CHECK_CLOSE(Acts::abs(distsAlongStrip[i] + truthDistances[i]), 0.,
                tolerance);
  }

  COMPARE_VECTORS(seedPos, truthPosZ0, tolerance);
  // check the direction
  CHECK_CLOSE(std::abs(seedDir.dot(truthDir)), 1.,
              std::numeric_limits<float>::epsilon());
}

void timeStripResidualTest(const Pars_t& linePars, const double timeT0,
                           const TestSpacePoint& sp,
                           const Acts::Transform3& locToGlob) {
  using ResidualIdx = CompSpacePointAuxiliaries::ResidualIdx;
  Config_t resCfg{};
  resCfg.localToGlobal = locToGlob;
  resCfg.useHessian = true;
  resCfg.calcAlongStrip = true;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::y0, ParIdx::phi, ParIdx::theta,
                      ParIdx::t0};
  Line_t line{linePars};

  ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Residual test for line: "
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

  CHECK_CLOSE(resCalc.residual()[toUnderlying(ResidualIdx::time)],
              (sp.time() - ToF), 1.e-10);
  ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Time of flight: " << ToF
                      << ", measured time : " << sp.time() << " --> residual: "
                      << (sp.time() - ToF) << " vs. from calculator: "
                      << resCalc.residual()[toUnderlying(ResidualIdx::time)]
                      << ".");
  constexpr double h = 5.e-9;

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
      lineUp.updateParameters(lineParsUp);
      lineDn.updateParameters(lineParsDn);
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
    COMPARE_VECTORS(numDeriv, resCalc.gradient(partial), tolerance);
    for (const auto partial2 : resCfg.parsToUse) {
      const Vector numDeriv1{
          (resCalcUp.gradient(partial2) - resCalcDn.gradient(partial2)) /
          (2. * h)};
      ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Second deriv ("
                          << CompSpacePointAuxiliaries::parName(partial) << ", "
                          << CompSpacePointAuxiliaries::parName(partial2)
                          << ") -- numerical: " << toString(numDeriv1)
                          << ", analytic: "
                          << toString(resCalc.hessian(partial, partial2)));
      COMPARE_VECTORS(numDeriv1, resCalc.hessian(partial, partial2), tolerance);
    }
  }
}
BOOST_AUTO_TEST_CASE(StrawDriftTimeCase) {
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
    ACTS_DEBUG(__func__ << "() " << __LINE__ << ": Calculate residual w.r.t. "
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
    ACTS_DEBUG(__func__ << "() " << __LINE__ << ": Create new space point @ "
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

      COMPARE_VECTORS(numDeriv, resCalc.gradient(partial), tolerance);
      for (const auto partial2 : resCfg.parsToUse) {
        const Vector numDeriv1{
            (resCalcUp.gradient(partial2) - resCalcDn.gradient(partial2)) /
            (2. * h)};
        ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Second deriv ("
                            << CompSpacePointAuxiliaries::parName(partial)
                            << ", "
                            << CompSpacePointAuxiliaries::parName(partial2)
                            << ") -- numerical: " << toString(numDeriv1)
                            << ", analytic: "
                            << toString(resCalc.hessian(partial, partial2)));
        COMPARE_VECTORS(numDeriv1, resCalc.hessian(partial, partial2),
                        tolerance);
      }
    }
  };

  RandomEngine rndEngine{4711};

  const Vector wirePos{100._cm, 50._cm, 30._cm};
  for (std::size_t e = 0; e < nEvents; ++e) {
    break;
    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Run test event: " << e);
    resCfg.localToGlobal = Acts::Transform3::Identity();
    Pars_t linePars{generateLine(rndEngine)};
    const double t0 = uniform{0_ns, 150._ns}(rndEngine);

    testTimingResidual(linePars, wirePos, Vector::UnitX(), t0);
    testTimingResidual(linePars, wirePos, Vector{1., 1., 0.}.normalized(), t0);

    resCfg.localToGlobal.translation() =
        Vector{uniform{-10._cm, 10._cm}(rndEngine),
               uniform{-20._cm, 20._cm}(rndEngine),
               uniform{-30._cm, 30._cm}(rndEngine)};
    testTimingResidual(linePars, wirePos, Vector::UnitX(), t0);
    testTimingResidual(linePars, wirePos, Vector{1., 1., 0.}.normalized(), t0);

    //// Next test the displacement
    resCfg.localToGlobal *=
        Acts::AngleAxis3{uniform{-30._degree, 30._degree}(rndEngine),
                         makeDirectionFromPhiTheta(
                             uniform{-30._degree, 30._degree}(rndEngine),
                             uniform{-45._degree, -35._degree}(rndEngine))};
    testTimingResidual(linePars, wirePos, Vector::UnitX(), t0);
    testTimingResidual(linePars, wirePos, Vector{1., 1., 0.}.normalized(), t0);
  }
}
BOOST_AUTO_TEST_CASE(WireResidualTest) {
  RandomEngine rndEngine{2525};
  ACTS_INFO("Run WireResidualTest");
  const Vector wirePos{100._cm, 50._cm, 30._cm};
  const Vector wireDir1{Vector::UnitX()};
  const Vector wireDir2{makeDirectionFromPhiTheta(10._degree, 30._degree)};

  for (std::size_t e = 0; e < nEvents; ++e) {
    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Run test event: " << e);
    Pars_t linePars{generateLine(rndEngine)};
    // Generate the first test measurement
    testResidual(linePars, TestSpacePoint{wirePos, wireDir1, 10._cm});

    testResidual(linePars, TestSpacePoint{wirePos, wireDir2, 10._cm});

    const Line_t line{linePars};
    auto isect =
        lineIntersect(line.position(), line.direction(), wirePos, wireDir1);
    testResidual(linePars, TestSpacePoint{
                               isect.position() +
                                   uniform{-50._cm, 50._cm}(rndEngine)*wireDir1,
                               wireDir1, 10._cm, true});
    isect = lineIntersect(line.position(), line.direction(), wirePos, wireDir2);
    testResidual(linePars, TestSpacePoint{
                               isect.position() +
                                   uniform{-50._cm, 50._cm}(rndEngine)*wireDir2,
                               wireDir2, 10._cm, true});
  }
}
BOOST_AUTO_TEST_CASE(StripResidual) {
  ACTS_INFO("Run StripResidualTest");
  RandomEngine rndEngine{2505};
  const Vector stripPos{75._cm, -75._cm, 100._cm};
  const Vector b1{makeDirectionFromPhiTheta(90._degree, 90._degree)};
  const Vector b2{makeDirectionFromPhiTheta(0._degree, 90_degree)};
  for (std::size_t e = 0; e < nEvents; ++e) {
    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Run test event: " << e);
    Pars_t linePars{generateLine(rndEngine)};

    testResidual(linePars,
                 TestSpacePoint{stripPos, b1, b2, TestSpacePoint::bendingDir});

    testResidual(linePars, TestSpacePoint{stripPos, b2, b1,
                                          TestSpacePoint::nonBendingDir});

    testResidual(linePars, TestSpacePoint{stripPos, b1, b2,
                                          TestSpacePoint::bothDirections});

    testResidual(
        linePars,
        TestSpacePoint{stripPos,
                       makeDirectionFromPhiTheta(30._degree, 90._degree),
                       makeDirectionFromPhiTheta(60._degree, 90_degree),
                       TestSpacePoint::bothDirections});
  }
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
    ACTS_DEBUG(__func__ << "() " << __LINE__ << ": Test space point "
                        << testMe);

    resCalc.updateFullResidual(line, t0, testMe);

    resCalc.updateChiSq(chi2, testMe.covariance());
    resCalc.symmetrizeHessian(chi2);
    constexpr auto N = 5u;
    for (std::size_t i = 1; i < N; ++i) {
      for (std::size_t j = 0; j <= i; ++j) {
        BOOST_CHECK_EQUAL(chi2.hessian(i, j), chi2.hessian(j, i));
      }
    }

    ACTS_DEBUG(__func__ << "() - " << __LINE__
                        << ": Calculated chi2: " << chi2.chi2 << ", Gradient: "
                        << toString(chi2.gradient) << ", Hessian: \n"
                        << chi2.hessian
                        << "\ndeterminant: " << chi2.hessian.determinant());
    CHECK_CLOSE(CompSpacePointAuxiliaries::chi2Term(line, resCfg.localToGlobal,
                                                    t0, testMe),
                chi2.chi2, 1.e-10);

    constexpr double h = 1.e-7;
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
      CHECK_CLOSE(numDeriv, anaDeriv, tolerance);
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
        CHECK_CLOSE(anaHess, numHess, tolerance);
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

  /// Test straws with information on the position along the
  /// tube
  testChi2(TestSpacePoint{line.point(20._cm) +
                              5._cm * line.direction().cross(Vector::UnitX()) +
                              4._cm * Vector::UnitX(),
                          Vector3::UnitX(),
                          5._cm,
                          true,
                          {Acts::pow(0.5_cm, 2), Acts::pow(10._mm, 2), 0.}});
}

BOOST_AUTO_TEST_CASE(CombinatorialSeedSolverStripsTest) {
  ACTS_INFO("Run Combinatorial seed test");
  RandomEngine rndEngine{23568};
  constexpr std::size_t nStrips = 8;

  const std::array<Vector3, nStrips> stripDirections{
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
    ACTS_DEBUG(__func__ << "() " << __LINE__ << " - Combinatorial Seed test - "
                        << "Processing Event: " << i);
    // update pseudo track parameters with random values
    Line_t line{generateLine(rndEngine)};

    const Vector3& muonPos = line.position();
    const Vector3& muonDir = line.direction();

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

}  // namespace ActsTests
BOOST_AUTO_TEST_SUITE_END()

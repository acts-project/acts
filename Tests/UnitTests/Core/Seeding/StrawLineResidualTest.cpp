// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/CombinatorialSeedSolver.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/detail/Polynomials.hpp"

#include "StrawHitGeneratorHelper.hpp"

using namespace Acts;
using namespace Acts::detail;
using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;
using namespace Acts::UnitLiterals;
using namespace Acts::PlanarHelper;
using namespace Acts::detail::LineHelper;
using namespace Acts::VectorHelpers;

using Line_t = CompSpacePointAuxiliaries::Line_t;
using Config_t = CompSpacePointAuxiliaries::Config;
using ParIdx = CompSpacePointAuxiliaries::FitParIndex;
using ResidualIdx = CompSpacePointAuxiliaries::ResidualIdx;

using Vector_t = Line_t::Vector;
using Pars_t = Line_t::ParamVector;

constexpr auto logLvl = Logging::Level::INFO;
constexpr std::size_t nEvents = 500;
constexpr double tolerance = 1.e-3;

ACTS_LOCAL_LOGGER(getDefaultLogger("StrawLineResidualTest", logLvl));

namespace ActsTests {

template <typename T>
constexpr T absMax(const T& a, const T& b) {
  return std::max(Acts::abs(a), Acts::abs(b));
}
template <typename T, typename... argsT>
constexpr T absMax(const T& a, const argsT... args) {
  return std::max(Acts::abs(a), absMax(args...));
}
double angle(const Vector_t& v1, const Vector_t& v2) {
  const double dots = Acts::abs(v1.dot(v2));
  return std::acos(std::clamp(dots, 0., 1.));
}

template <int N>
bool isPositiveDefinite(const SquareMatrix<N>& mat) {
  if (mat.determinant() < Acts::s_epsilon) {
    return false;
  }
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j <= i; ++j) {
      if (mat(i, j) != mat(j, i)) {
        return false;
      }
    }
  }
  Eigen::EigenSolver<Acts::SquareMatrix<N>> eigenDecomp{mat};
  for (const auto lambda : eigenDecomp.eigenvalues()) {
    if (Acts::abs(lambda.imag()) > Acts::s_epsilon ||
        lambda.real() < Acts::s_epsilon) {
      return false;
    }
  }
  return true;
}

std::string whiteSpaces(const std::size_t n, const char s = '#') {
  std::string str{};
  str.resize(n, s);
  return str;
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

void testResidual(const Pars_t& linePars, const FitTestSpacePoint& testPoint) {
  Config_t resCfg{};
  resCfg.useHessian = true;
  resCfg.calcAlongStrip = false;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::y0, ParIdx::phi, ParIdx::theta};
  Line_t line{linePars};

  ACTS_INFO(__func__ << "() - " << __LINE__ << ":\n"
                     << whiteSpaces(150) << "\n"
                     << whiteSpaces(150) << "\nResidual test w.r.t. test line: "
                     << toString(line.position()) << ", "
                     << toString(line.direction())
                     << ", theta: " << (theta(line.direction()) / 1._degree)
                     << ", phi: " << (phi(line.direction())) / 1._degree);

  CompSpacePointAuxiliaries resCalc{resCfg,
                                    Acts::getDefaultLogger("testRes", logLvl)};
  CompSpacePointAuxiliaries resCalcUp{
      resCfg, Acts::getDefaultLogger("testResUp", logLvl)};
  CompSpacePointAuxiliaries resCalcDn{
      resCfg, Acts::getDefaultLogger("testResDn", logLvl)};

  resCalc.updateSpatialResidual(line, testPoint);
  ACTS_INFO(__func__ << "() - " << __LINE__ << ": Test residual w.r.t. "
                     << toString(testPoint));
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
                  closePoint.pathLength(), 1.e-9);
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

    Acts::SquareMatrix<3> coordTrf{Acts::SquareMatrix<3>::Zero()};
    coordTrf.col(toUnderlying(ResidualIdx::bending)) =
        testPoint.measuresLoc1() ? testPoint.toNextSensor()
                                 : testPoint.sensorDirection();
    coordTrf.col(toUnderlying(ResidualIdx::nonBending)) =
        testPoint.measuresLoc1() ? testPoint.sensorDirection()
                                 : testPoint.toNextSensor();
    coordTrf.col(2) = n;
    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Coordinate trf: \n"
                        << toString(coordTrf));

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

    resCalcUp.reset();
    resCalcDn.reset();

    resCalcUp.updateSpatialResidual(lineUp, testPoint);
    resCalcDn.updateSpatialResidual(lineDn, testPoint);

    const Vector_t numDeriv =
        (resCalcUp.residual() - resCalcDn.residual()) / (2. * h);

    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Derivative test: "
                        << CompSpacePointAuxiliaries::parName(par)
                        << ", derivative: " << toString(resCalc.gradient(par))
                        << ",  numerical: " << toString(numDeriv));

    COMPARE_VECTORS(numDeriv, resCalc.gradient(par), tolerance);
    /// Next attempt to calculate the second derivative
    for (auto par1 : resCfg.parsToUse) {
      const Vector_t numDeriv1{
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

template <typename Point_t, std::size_t N>
void testSeed(const std::array<Point_t, N>& spacePoints,
              const std::array<double, N>& truthDistances,
              const Vector3& truthPosZ0, const Vector3& truthDir) {
  const SquareMatrix2 bMatrix =
      CombinatorialSeedSolver::betaMatrix(spacePoints);

  if (std::abs(bMatrix.determinant()) < std::numeric_limits<float>::epsilon()) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": Beta Matrix \n"
                          << toString(bMatrix)
                          << "\nhas zero determinant - skip the combination");
    return;
  }

  const auto distsAlongStrip =
      CombinatorialSeedSolver::defineParameters(bMatrix, spacePoints);

  const auto [seedPos, seedDir] =
      CombinatorialSeedSolver::seedSolution(spacePoints, distsAlongStrip);

  ACTS_DEBUG(__func__ << "() " << __LINE__ << ": Reconstructed position "
                      << toString(seedPos) << " and direction "
                      << toString(seedDir));
  ACTS_DEBUG(__func__ << "() " << __LINE__ << ": Truth position "
                      << toString(truthPosZ0) << " and truth direction "
                      << toString(truthDir));

  for (std::size_t i = 0; i < distsAlongStrip.size(); ++i) {
    CHECK_CLOSE(Acts::abs(distsAlongStrip[i] + truthDistances[i]), 0.,
                tolerance);
  }

  COMPARE_VECTORS(seedPos, truthPosZ0, tolerance);
  CHECK_CLOSE(std::abs(seedDir.dot(truthDir)), 1.,
              std::numeric_limits<float>::epsilon());
}

void timeStripResidualTest(const Pars_t& linePars, const double timeT0,
                           const FitTestSpacePoint& sp,
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

    const Vector_t numDeriv =
        (resCalcUp.residual() - resCalcDn.residual()) / (2. * h);

    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Derivative test: "
                        << CompSpacePointAuxiliaries::parName(partial)
                        << ", derivative: "
                        << toString(resCalc.gradient(partial))
                        << ",  numerical: " << toString(numDeriv));
    COMPARE_VECTORS(numDeriv, resCalc.gradient(partial), tolerance);
    for (const auto partial2 : resCfg.parsToUse) {
      const Vector_t numDeriv1{
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

BOOST_AUTO_TEST_CASE(CalibratorTest) {
  constexpr double h = 1.e-7 * 1._ns;
  for (double t = detailCalib::s_minDriftTime + 1._ns;
       t < detailCalib::s_maxDriftTime; t += 10._ns) {
    const double r = detailCalib::driftRadius(t);
    const double v = detailCalib::driftVelocity(t);
    const double a = detailCalib::driftAcceleration(t);

    ACTS_DEBUG(__func__ << "() " << __LINE__ << " - t:" << (t / 1._ns)
                        << " r: " << r << ", v: " << v << ", a: " << a);

    const double numV =
        (detailCalib::driftRadius(t + h) - detailCalib::driftRadius(t - h)) /
        (2. * h);
    const double numA = (detailCalib::driftVelocity(t + h) -
                         detailCalib::driftVelocity(t - h)) /
                        (2. * h);

    CHECK_CLOSE(numV, v, tolerance);
    CHECK_CLOSE(numA, a, tolerance);

    ACTS_DEBUG(__func__ << "() " << __LINE__ << ": numV: " << numV
                        << ", numA: " << numA);

    const double backT = detailCalib::driftTime(r);
    const double backTprime = detailCalib::driftTimePrime(r);

    const double numTprime =
        (detailCalib::driftTime(r + h) - detailCalib::driftTime(r - h)) /
        (2. * h);
    CHECK_CLOSE(numTprime, backTprime, tolerance);
    ACTS_DEBUG(__func__ << "() " << __LINE__ << " - backT: " << (backT / 1._ns)
                        << ", back T^{prime}: " << (backTprime / 1._ns));
  }
}

BOOST_AUTO_TEST_CASE(StrawsWithT0Constraint) {
  ACTS_INFO("Calibration straw point test ");

  Config_t resCfg{};
  resCfg.useHessian = true;
  resCfg.calcAlongStrip = false;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::phi, ParIdx::y0, ParIdx::theta,
                      ParIdx::t0};

  CompSpacePointAuxiliaries resCalc{resCfg,
                                    Acts::getDefaultLogger("testRes", logLvl)};
  CompSpacePointAuxiliaries resCalcUp{
      resCfg, Acts::getDefaultLogger("testResUp", logLvl)};
  CompSpacePointAuxiliaries resCalcDn{
      resCfg, Acts::getDefaultLogger("testResDn", logLvl)};

  MeasurementGenerator::Config genCfg{};
  genCfg.createStraws = true;
  genCfg.smearRadius = false;
  genCfg.twinStraw = false;
  genCfg.createStrips = false;

  RandomEngine rndEngine{4711};

  constexpr double h = 5.e-8;

  const Acts::CalibrationContext cctx{};
  SpCalibrator calibrator{logLvl};

  for (std::size_t e = 0; e < nEvents; ++e) {
    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Run test event: " << e);
    resCfg.localToGlobal = Acts::Transform3::Identity();
    Line_t line{generateLine(rndEngine, logger())};
    const Pars_t linePars = line.parameters();
    const double t0 = uniform{0_ns, 50._ns}(rndEngine);

    const auto container =
        MeasurementGenerator::spawn(line, t0, rndEngine, genCfg, logger());
    BOOST_CHECK_GE(container.size(), 1u);
    const auto testMeas = container.front();

    resCalc.reset();

    resCalc.updateFullResidual(line, t0, *testMeas,
                               calibrator.driftVelocity(cctx, *testMeas),
                               calibrator.driftAcceleration(cctx, *testMeas));

    for (const auto par : resCfg.parsToUse) {
      Pars_t lineParsUp{linePars};
      Pars_t lineParsDn{linePars};

      const auto dIdx = toUnderlying(par);
      double t0Up{t0}, t0Dn{t0};
      if (par != ParIdx::t0) {
        lineParsUp[dIdx] += h;
        lineParsDn[dIdx] -= h;
      } else {
        t0Up += h;
        t0Dn -= h;
      }

      resCalcUp.reset();
      line.updateParameters(lineParsUp);
      const auto testMeasUp = calibrator.calibrate(
          cctx, line.position(), line.direction(), t0Up, *testMeas);
      resCalcUp.updateFullResidual(
          line, t0Up, *testMeasUp, calibrator.driftVelocity(cctx, *testMeasUp),
          calibrator.driftAcceleration(cctx, *testMeasUp));

      resCalcDn.reset();
      line.updateParameters(lineParsDn);
      const auto testMeasDn = calibrator.calibrate(
          cctx, line.position(), line.direction(), t0Dn, *testMeas);

      resCalcDn.updateFullResidual(
          line, t0Dn, *testMeasDn, calibrator.driftVelocity(cctx, *testMeasDn),
          calibrator.driftAcceleration(cctx, *testMeasDn));

      const Vector_t numDeriv =
          (resCalcUp.residual() - resCalcDn.residual()) / (2. * h);
      const Vector_t anaDeriv = resCalc.gradient(par);
      COMPARE_VECTORS(numDeriv, anaDeriv, 1.e-3);
      for (const auto par1 : resCfg.parsToUse) {
        const Vector_t numDeriv1 =
            (resCalcUp.gradient(par1) - resCalcDn.gradient(par1)) / (2. * h);
        const Vector_t anaDeriv1 = resCalc.hessian(par1, par);
        COMPARE_VECTORS(numDeriv1, anaDeriv1, 4.e-3);
      }
    }
  }
}
BOOST_AUTO_TEST_CASE(WireResidualTest) {
  RandomEngine rndEngine{2525};
  ACTS_INFO("Run WireResidualTest");
  const Vector_t wirePos{100._cm, 50._cm, 30._cm};
  const Vector_t wireDir1{Vector_t::UnitX()};
  const Vector_t wireDir2{makeDirectionFromPhiTheta(10._degree, 30._degree)};

  for (std::size_t e = 0; e < nEvents; ++e) {
    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Run test event: " << e);
    const Line_t line{generateLine(rndEngine, logger())};
    // Generate the first test measurement
    constexpr double R = 13._mm;
    const double dR = detailCalib::driftRadUncert(R);
    constexpr double uncertWire = 1._cm;
    testResidual(line.parameters(),
                 FitTestSpacePoint{wirePos, wireDir1, R, dR});

    testResidual(line.parameters(),
                 FitTestSpacePoint{wirePos, wireDir2, R, dR});

    auto isect =
        lineIntersect(line.position(), line.direction(), wirePos, wireDir1);

    testResidual(line.parameters(),
                 FitTestSpacePoint{isect.position() + uniform{-50._cm, 50._cm}(
                                                          rndEngine)*wireDir1,
                                   wireDir1, R, dR, uncertWire});

    isect = lineIntersect(line.position(), line.direction(), wirePos, wireDir2);
    testResidual(line.parameters(),
                 FitTestSpacePoint{isect.position() + uniform{-50._cm, 50._cm}(
                                                          rndEngine)*wireDir2,
                                   wireDir2, R, dR, uncertWire});
  }
}
BOOST_AUTO_TEST_CASE(StripResidual) {
  ACTS_INFO("Run StripResidualTest");
  RandomEngine rndEngine{2505};
  const Vector_t stripPos{75._cm, -75._cm, 100._cm};
  const Vector_t b1{makeDirectionFromPhiTheta(90._degree, 90._degree)};
  const Vector_t b2{makeDirectionFromPhiTheta(0._degree, 90_degree)};
  for (std::size_t e = 0; e < nEvents; ++e) {
    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Run test event: " << e);
    Pars_t linePars{generateLine(rndEngine, logger()).parameters()};

    /// Test solely the residual along the bending direction
    testResidual(linePars, FitTestSpacePoint{stripPos, b1, b2, 0._cm, 1._cm});
    /// Test solely the residual along the non-bending direction
    testResidual(linePars, FitTestSpacePoint{stripPos, b2, b1, 1._cm, 0._cm});
    /// Test the combined residual
    testResidual(linePars, FitTestSpacePoint{stripPos, b1, b2, 1._cm, 1._cm});
    const Vector_t b3 = makeDirectionFromPhiTheta(30._degree, 90._degree);
    const Vector_t b4 = makeDirectionFromPhiTheta(60._degree, 90._degree);

    testResidual(linePars, FitTestSpacePoint{stripPos, b3, b4, 1._cm, 1._cm});
  }
}

BOOST_AUTO_TEST_CASE(TimeStripResidual) {
  Pars_t linePars{};
  linePars[toUnderlying(ParIdx::phi)] = 60._degree;
  linePars[toUnderlying(ParIdx::theta)] = 45_degree;

  Acts::Transform3 locToGlob{Acts::Transform3::Identity()};

  const Vector_t pos{75._cm, -75._cm, 100._cm};
  const Vector_t pos1{75._cm, -75._cm, 200._cm};

  const Vector_t b1{makeDirectionFromPhiTheta(30._degree, 90._degree)};
  const Vector_t b2{makeDirectionFromPhiTheta(60._degree, 90._degree)};
  const Vector_t b3{makeDirectionFromPhiTheta(00._degree, 90._degree)};
  const Vector_t b4{makeDirectionFromPhiTheta(60._degree, 75._degree)};
  const std::array cov{Acts::square(10._cm), Acts::square(10._cm),
                       Acts::square(1._ns)};
  FitTestSpacePoint p1{pos, b1, b2, 15._ns, cov};

  FitTestSpacePoint p2{pos1, b3, b4, 15._ns, cov};

  timeStripResidualTest(linePars, 10., p1, locToGlob);

  timeStripResidualTest(linePars, 10., p2, locToGlob);

  locToGlob.translation() = Vector_t{75._cm, -75._cm, -35._cm};
  timeStripResidualTest(linePars, 10., p1, locToGlob);

  timeStripResidualTest(linePars, 10., p2, locToGlob);
  locToGlob *= Acts::AngleAxis3{
      30_degree, makeDirectionFromPhiTheta(30_degree, -45_degree)};
  timeStripResidualTest(linePars, 10., p1, locToGlob);

  timeStripResidualTest(linePars, 10., p2, locToGlob);
}

BOOST_AUTO_TEST_CASE(ChiSqEvaluation) {
  Config_t resCfg{};
  resCfg.useHessian = true;
  resCfg.calcAlongStrip = false;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::phi, ParIdx::y0, ParIdx::theta,
                      ParIdx::t0};

  using ResCalc_t = CompSpacePointAuxiliaries;
  using ChiSq_t = ResCalc_t::ChiSqWithDerivatives;

  ResCalc_t resCalc{resCfg, Acts::getDefaultLogger("chiSquareTest", logLvl)};

  ChiSq_t chi2{};
  ChiSq_t chi2Up{};
  ChiSq_t chi2Dn{};

  Line_t line{};

  const Acts::CalibrationContext cctx{};

  const SpCalibrator calibrator{logLvl};

  RandomEngine rndEngine{1789};

  auto runCalculation = [&line, &resCalc, &calibrator, &cctx](
                            FitTestSpacePoint testMe, const Pars_t& pars,
                            const double t0, ChiSq_t& chiSq, const bool calib) {
    line.updateParameters(pars);
    resCalc.reset();
    chiSq.reset();
    if (calib && testMe.isStraw() &&
        Acts::rangeContainsValue(resCalc.config().parsToUse, ParIdx::t0)) {
      testMe = *calibrator.calibrate(cctx, line.position(), line.direction(),
                                     t0, testMe);
    }

    if (testMe.isStraw()) {
      const int sign = CompSpacePointAuxiliaries::strawSign(line, testMe);
      testMe.updateDriftR(std::copysign(testMe.driftRadius(), sign));
    }

    if (Acts::rangeContainsValue(resCalc.config().parsToUse, ParIdx::t0)) {
      resCalc.updateFullResidual(line, t0, testMe,
                                 calibrator.driftVelocity(cctx, testMe),
                                 calibrator.driftAcceleration(cctx, testMe));
    } else {
      resCalc.updateSpatialResidual(line, testMe);
    }
    resCalc.updateChiSq(chiSq, testMe.covariance());
    resCalc.symmetrizeHessian(chiSq);
  };

  auto testChi2 = [&](const Pars_t& linePars, const double t0,
                      FitTestSpacePoint testMe) {
    ACTS_DEBUG(__func__ << "() " << __LINE__ << ":\n"
                        << whiteSpaces(150, '+') << "\nTest space point "
                        << toString(testMe) << "\n"
                        << whiteSpaces(150, '+'));
    line.updateParameters(linePars);
    double chi2Term{0.};
    if (Acts::rangeContainsValue(resCfg.parsToUse, ParIdx::t0)) {
      chi2Term = CompSpacePointAuxiliaries::chi2Term(line, resCfg.localToGlobal,
                                                     t0, testMe);
    } else {
      chi2Term = CompSpacePointAuxiliaries::chi2Term(line, testMe);
    }
    if (testMe.isStraw()) {
      testMe.updateTime(detailCalib::driftTime(testMe.driftRadius()) + t0 +
                        calibrator.fastToF(testMe));
    }

    ACTS_DEBUG(__func__ << "() - " << __LINE__
                        << " Evaluated chi2 term: " << chi2Term);
    runCalculation(testMe, linePars, t0, chi2, false);
    constexpr auto N = 5u;
    for (std::size_t i = 1; i < N; ++i) {
      for (std::size_t j = 0; j <= i; ++j) {
        BOOST_CHECK_EQUAL(chi2.hessian(i, j), chi2.hessian(j, i));
      }
    }

    ACTS_DEBUG(__func__ << "() - " << __LINE__
                        << ": Calculated chi2: " << chi2.chi2 << ", Gradient: "
                        << toString(chi2.gradient) << ", Hessian: \n"
                        << toString(chi2.hessian)
                        << "\ndeterminant: " << chi2.hessian.determinant());
    CHECK_CLOSE(chi2Term, chi2.chi2, 1.e-10);
    runCalculation(testMe, linePars, t0, chi2, true);
    constexpr double h = 5.e-8;
    for (const auto par : resCfg.parsToUse) {
      Pars_t lineParsUp{linePars};
      Pars_t lineParsDn{linePars};

      const auto dIdx = toUnderlying(par);
      double t0Up{t0}, t0Dn{t0};
      if (par != ParIdx::t0) {
        lineParsUp[dIdx] += h;
        lineParsDn[dIdx] -= h;
      } else {
        t0Up += h;
        t0Dn -= h;
      }
      runCalculation(testMe, lineParsUp, t0Up, chi2Up, true);
      runCalculation(testMe, lineParsDn, t0Dn, chi2Dn, true);

      const double anaDeriv = chi2.gradient[dIdx];
      const double numDeriv = (chi2Up.chi2 - chi2Dn.chi2) / (2. * h);

      ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Partial "
                          << CompSpacePointAuxiliaries::parName(par)
                          << " --  numerical: " << numDeriv
                          << ", analytical: " << anaDeriv);
      CHECK_CLOSE(numDeriv, anaDeriv, tolerance);
      for (const auto par2 : resCfg.parsToUse) {
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
    ACTS_DEBUG(__func__ << "() " << __LINE__ << ":\n"
                        << whiteSpaces(150, '-') << "\n"
                        << whiteSpaces(150, '-'));
  };

  const double t0{2.5_ns};

  for (std::size_t e = 0; e < nEvents; ++e) {
    line = generateLine(rndEngine, logger());
    const Pars_t pars = line.parameters();

    /// Test orthogonal strips
    testChi2(
        pars, t0,
        FitTestSpacePoint{
            line.point(20._cm) +
                5._cm * line.direction().cross(Vector_t::UnitX()),
            makeDirectionFromPhiTheta(0_degree, 90._degree),
            makeDirectionFromPhiTheta(90._degree, 0._degree),
            15._ns,
            {Acts::pow(5._cm, 2), Acts::pow(10._cm, 2), Acts::pow(1._ns, 2)}});

    // Test strips with stereo
    testChi2(
        pars, t0,
        FitTestSpacePoint{
            line.point(20._cm) +
                5._cm * line.direction().cross(Vector_t::UnitX()),
            makeDirectionFromPhiTheta(0_degree, 45._degree),
            makeDirectionFromPhiTheta(60._degree, 0._degree),
            15._ns,
            {Acts::pow(5._cm, 2), Acts::pow(10._cm, 2), Acts::pow(1._ns, 2)}});
    constexpr double R = 12._mm;
    const double dR = detailCalib::driftRadUncert(R);
    // Test ordinary straws
    testChi2(pars, t0,
             FitTestSpacePoint{line.point(20._cm) + R * line.direction().cross(
                                                            Vector_t::UnitX()),
                               Vector3::UnitX(), R, dR});

    testChi2(pars, t0,
             FitTestSpacePoint{line.point(40._cm) - R * line.direction().cross(
                                                            Vector_t::UnitX()),
                               Vector3::UnitX(), R, dR});

    // Test straws with information on the position along the tube
    constexpr double dZ = 2._cm;
    testChi2(
        pars, t0,
        FitTestSpacePoint{line.point(10._cm) +
                              R * line.direction().cross(Vector_t::UnitX()) +
                              dZ * Vector_t::UnitX(),
                          Vector3::UnitX(), R, dR, dZ});
    testChi2(
        pars, t0,
        FitTestSpacePoint{line.point(10._cm) -
                              R * line.direction().cross(Vector_t::UnitX()) +
                              dZ * Vector_t::UnitX(),
                          Vector3::UnitX(), R, dR, dZ});
  }
}

BOOST_AUTO_TEST_CASE(HessianCalc) {
  // Check that the hessian matrices are positive definite
  Config_t resCfg{};
  resCfg.useHessian = false;
  resCfg.calcAlongStrip = true;
  resCfg.parsToUse = {ParIdx::y0, ParIdx::theta};

  CompSpacePointAuxiliaries resCalc{
      resCfg, Acts::getDefaultLogger("ResidualCalcTwin", logLvl)};
  RandomEngine rndEngine{1167};

  using ChiSq_t = CompSpacePointAuxiliaries::ChiSqWithDerivatives;
  ChiSq_t chi2{};

  MeasurementGenerator::Config genCfg{};
  genCfg.createStraws = true;
  genCfg.smearRadius = true;
  genCfg.twinStraw = false;
  genCfg.createStrips = false;

  for (std::size_t ev = 0; ev < nEvents; ++ev) {
    break;
    resCalc.reset();
    chi2.reset();
    const Line_t line = generateLine(rndEngine, logger());
    for (const auto& tube :
         MeasurementGenerator::spawn(line, 0., rndEngine, genCfg, logger())) {
      resCalc.updateSpatialResidual(line, *tube);
      resCalc.updateChiSq(chi2, tube->covariance());
    }
    resCalc.symmetrizeHessian(chi2);

    const SquareMatrix<2> miniHessian{chi2.hessian.block<2, 2>(0, 0)};
    BOOST_CHECK_EQUAL(isPositiveDefinite(miniHessian), true);

    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ":\n"
                        << miniHessian << ", determinant: "
                        << miniHessian.determinant() << ",\n"
                        << printEigenDecomposition(miniHessian));
  }
}

BOOST_AUTO_TEST_CASE(HessianCalcT0) {
  // Check that the hessian matrices are positive definite
  Config_t resCfg{};
  resCfg.useHessian = false;
  resCfg.calcAlongStrip = true;
  resCfg.parsToUse = {ParIdx::y0, ParIdx::theta, ParIdx::t0};

  CompSpacePointAuxiliaries resCalc{
      resCfg, Acts::getDefaultLogger("ResidualCalcTwin", logLvl)};
  RandomEngine rndEngine{1167};

  using ChiSq_t = CompSpacePointAuxiliaries::ChiSqWithDerivatives;
  ChiSq_t chi2{};

  MeasurementGenerator::Config genCfg{};
  genCfg.createStraws = true;
  genCfg.smearRadius = true;
  genCfg.twinStraw = false;
  genCfg.createStrips = false;

  SpCalibrator calibrator{logLvl};
  CalibrationContext cctx{};

  for (std::size_t ev = 0; ev < nEvents; ++ev) {
    resCalc.reset();
    chi2.reset();
    const Line_t line = generateLine(rndEngine, logger());
    const double t0 = 0.;
    for (const auto& tube :
         MeasurementGenerator::spawn(line, t0, rndEngine, genCfg, logger())) {
      resCalc.updateFullResidual(line, t0, *tube,
                                 calibrator.driftVelocity(cctx, *tube));
      resCalc.updateChiSq(chi2, tube->covariance());
    }
    resCalc.symmetrizeHessian(chi2);

    chi2.hessian.col(2).swap(chi2.hessian.col(toUnderlying(ParIdx::t0)));
    chi2.hessian.row(2).swap(chi2.hessian.row(toUnderlying(ParIdx::t0)));

    const SquareMatrix<3> miniHessian{chi2.hessian.block<3, 3>(0, 0)};
    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ":\n"
                        << miniHessian << ", determinant: "
                        << miniHessian.determinant() << ",\n"
                        << printEigenDecomposition(miniHessian));
    BOOST_CHECK_EQUAL(isPositiveDefinite(miniHessian), true);
    assert(isPositiveDefinite(miniHessian));
  }
}

BOOST_AUTO_TEST_CASE(HessianCalcTwin) {
  // Check that the hessian matrices are positive definite
  Config_t resCfg{};
  resCfg.useHessian = false;
  resCfg.calcAlongStrip = true;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::phi, ParIdx::y0, ParIdx::theta};

  CompSpacePointAuxiliaries resCalc{
      resCfg, Acts::getDefaultLogger("ResidualCalcTwin", logLvl)};
  RandomEngine rndEngine{1066};

  using ChiSq_t = CompSpacePointAuxiliaries::ChiSqWithDerivatives;
  ChiSq_t chi2{};

  MeasurementGenerator::Config genCfg{};
  genCfg.createStraws = true;
  genCfg.smearRadius = true;
  genCfg.twinStraw = true;
  genCfg.smearTwin = true;
  genCfg.createStrips = false;

  for (std::size_t ev = 0; ev < nEvents; ++ev) {
    resCalc.reset();
    chi2.reset();
    const Line_t line = generateLine(rndEngine, logger());
    for (const auto& twinTube :
         MeasurementGenerator::spawn(line, 0., rndEngine, genCfg, logger())) {
      resCalc.updateSpatialResidual(line, *twinTube);
      resCalc.updateChiSq(chi2, twinTube->covariance());
    }
    resCalc.symmetrizeHessian(chi2);

    const SquareMatrix<4> miniHessian{chi2.hessian.block<4, 4>(0, 0)};
    BOOST_CHECK_EQUAL(isPositiveDefinite(miniHessian), true);

    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ":\n"
                        << miniHessian << ", determinant: "
                        << miniHessian.determinant() << ",\n"
                        << printEigenDecomposition(miniHessian));
  }
}

BOOST_AUTO_TEST_CASE(HessianCalcStrawStrip) {
  // Check that the hessian matrices are positive definite
  Config_t resCfg{};
  resCfg.useHessian = false;
  resCfg.calcAlongStrip = true;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::phi, ParIdx::y0, ParIdx::theta};

  CompSpacePointAuxiliaries resCalc{
      resCfg, Acts::getDefaultLogger("ResidualCalcStrawStrip", logLvl)};
  RandomEngine rndEngine{1905};

  using ChiSq_t = CompSpacePointAuxiliaries::ChiSqWithDerivatives;
  ChiSq_t chi2{};

  MeasurementGenerator::Config genCfg{};
  genCfg.createStraws = true;
  genCfg.smearRadius = true;
  genCfg.twinStraw = false;
  genCfg.createStrips = true;

  for (std::size_t ev = 0; ev < nEvents; ++ev) {
    resCalc.reset();
    chi2.reset();
    const Line_t line = generateLine(rndEngine, logger());
    for (const auto& meas :
         MeasurementGenerator::spawn(line, 0., rndEngine, genCfg, logger())) {
      resCalc.updateSpatialResidual(line, *meas);
      resCalc.updateChiSq(chi2, meas->covariance());
    }
    resCalc.symmetrizeHessian(chi2);

    const SquareMatrix<4> miniHessian{chi2.hessian.block<4, 4>(0, 0)};
    BOOST_CHECK_EQUAL(isPositiveDefinite(miniHessian), true);

    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ":\n"
                        << miniHessian << ", determinant: "
                        << miniHessian.determinant() << ",\n"
                        << printEigenDecomposition(miniHessian));
  }
}

BOOST_AUTO_TEST_CASE(HessianCalcStripOnly) {
  // Check that the hessian matrices are positive definite
  Config_t resCfg{};
  resCfg.useHessian = false;
  resCfg.calcAlongStrip = true;
  resCfg.parsToUse = {ParIdx::x0, ParIdx::phi, ParIdx::y0, ParIdx::theta};

  CompSpacePointAuxiliaries resCalc{
      resCfg, Acts::getDefaultLogger("ResidualCalcStrip", logLvl)};
  RandomEngine rndEngine{1492};

  using ChiSq_t = CompSpacePointAuxiliaries::ChiSqWithDerivatives;
  ChiSq_t chi2{};

  MeasurementGenerator::Config genCfg{};
  genCfg.createStraws = false;
  genCfg.smearRadius = true;
  genCfg.twinStraw = false;
  genCfg.createStrips = true;

  for (std::size_t ev = 0; ev < nEvents; ++ev) {
    resCalc.reset();
    chi2.reset();
    const Line_t line = generateLine(rndEngine, logger());
    for (const auto& strip :
         MeasurementGenerator::spawn(line, 0., rndEngine, genCfg, logger())) {
      resCalc.updateSpatialResidual(line, *strip);
      resCalc.updateChiSq(chi2, strip->covariance());
    }
    resCalc.symmetrizeHessian(chi2);

    const SquareMatrix<4> miniHessian{chi2.hessian.block<4, 4>(0, 0)};
    BOOST_CHECK_EQUAL(isPositiveDefinite(miniHessian), true);

    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ":\n"
                        << miniHessian << ", determinant: "
                        << miniHessian.determinant() << ",\n"
                        << printEigenDecomposition(miniHessian));
  }
}

BOOST_AUTO_TEST_CASE(CombinatorialSeedSolverStripsTest) {
  ACTS_INFO("Run Combinatorial seed test");
  RandomEngine rndEngine{23568};
  constexpr std::size_t nStrips = 8;
  constexpr std::size_t nCombined = 2;

  const std::array<Vector3, nStrips> stripDirections{
      Vector3::UnitX(),
      Vector3::UnitX(),
      makeDirectionFromPhiTheta(1.5_degree, 90_degree),
      makeDirectionFromPhiTheta(-1.5_degree, 90_degree),
      makeDirectionFromPhiTheta(1.5_degree, 90_degree),
      makeDirectionFromPhiTheta(-1.5_degree, 90_degree),
      Vector3::UnitX(),
      Vector3::UnitX()};

  // place the strip measuremens and the 2D combined measurements - one layer
  // before and one layer after the strips
  constexpr std::array<double, nStrips> distancesZ{-20., -264., 0.,  300.,
                                                   350,  370.,  400, 450.};

  constexpr std::array<double, nCombined> distancesZCombined{-350., 360.};

  std::array<double, nStrips> parameters{};
  std::array<Vector3, nStrips> intersections{};

  // pseudo track initialization

  for (std::size_t i = 0; i < nEvents; ++i) {
    ACTS_DEBUG(__func__ << "() " << __LINE__ << " - Combinatorial Seed test - "
                        << "Processing Event: " << i);
    // update pseudo track parameters with random values
    const Line_t line{generateLine(rndEngine, logger())};

    const Vector3& muonPos = line.position();
    const Vector3& muonDir = line.direction();

    ACTS_DEBUG("Seed position and direction: " << toString(muonPos) << ", "
                                               << toString(muonDir));

    std::array<std::shared_ptr<FitTestSpacePoint>, nStrips> spacePoints{};
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

      spacePoints[layer] = std::make_shared<FitTestSpacePoint>(
          stripPos, stripDirections[layer], Vector3::UnitZ(), 1._cm, 0._cm);
    }

    // construct the combined space points
    std::array<std::shared_ptr<FitTestSpacePoint>, nCombined> combinedPoints{};
    for (std::size_t c = 0; c < nCombined; ++c) {
      auto intersection = PlanarHelper::intersectPlane(
          muonPos, muonDir, Vector3::UnitZ(), distancesZCombined[c]);

      combinedPoints[c] = std::make_shared<FitTestSpacePoint>(
          intersection.position(), Vector3::Zero(), Vector3::UnitZ(), 1._cm,
          1._cm);
    }

    // let's pick four strips on four different layers
    // check combinatorics
    std::array<std::shared_ptr<FitTestSpacePoint>, 4> seedSpacePoints{};
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

    // test the three layers combinatorics
    std::array<std::shared_ptr<FitTestSpacePoint>, 3> seedTriplet{};
    for (unsigned int l = 0; l < distancesZ.size() - 1; ++l) {
      seedTriplet[0] = spacePoints[l];
      for (unsigned int k = l + 1; k < distancesZ.size(); ++k) {
        seedTriplet[1] = spacePoints[k];
        for (unsigned int c = 0; c < nCombined; ++c) {
          seedTriplet[2] = combinedPoints[c];

          const std::array<double, 3> truthDistances = {
              parameters[l], parameters[k],
              0.};  // assume zero distance for the combined spacepoint

          testSeed(seedTriplet, truthDistances, muonPos, muonDir);
        }
      }
    }
  }
}

}  // namespace ActsTests
BOOST_AUTO_TEST_SUITE_END()

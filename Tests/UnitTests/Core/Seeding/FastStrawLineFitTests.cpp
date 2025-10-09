// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/detail/FastStrawLineFitter.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <format>
#include <random>

#include "TFile.h"
#include "TTree.h"

using namespace Acts;
using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;
using namespace Acts::UnitLiterals;
using RandomEngine = std::mt19937;
using uniform_t = std::uniform_real_distribution<double>;
using gauss_t = std::normal_distribution<double>;

constexpr std::size_t nTrials = 1;
constexpr auto logLvl = Logging::Level::INFO;

namespace ActsTests {

ACTS_LOCAL_LOGGER(getDefaultLogger("FastStrawLineFitTests", logLvl));

class StrawTestPoint;
using TestStrawCont_t = std::vector<std::unique_ptr<StrawTestPoint>>;
using Line_t = CompSpacePointAuxiliaries::Line_t;
using ResidualIdx = FastStrawLineFitter::ResidualIdx;

/// @brief Outstream operator for a vector
template <typename T>
std::ostream& operator<<(std::ostream& ostr, const std::vector<T>& v) {
  ostr << "[";
  for (std::size_t i = 0; i < v.size(); ++i) {
    ostr << v[i];
    if (i + 1 != v.size()) {
      ostr << ", ";
    }
  }
  ostr << "]";
  return ostr;
}
/// @brief Converts an ACTS time into SI [ns]
constexpr double inNanoS(const double x) {
  return x / 1._ns;
}

class StrawTestPoint {
 public:
  StrawTestPoint(const Vector3& pos, const double driftR,
                 const double driftRUncert)
      : m_pos{pos}, m_driftR{abs(driftR)} {
    m_cov[toUnderlying(ResidualIdx::bending)] = square(driftRUncert);
  }

  StrawTestPoint(const Vector3& pos, const double uncert)
      : m_pos{pos}, m_isStraw{false} {
    m_cov[toUnderlying(ResidualIdx::bending)] = square(uncert);
  }
  /// @brief Straw tube's direction
  const Vector3& localPosition() const { return m_pos; }
  /// @brief Wire direction
  const Vector3& sensorDirection() const { return m_wireDir; }
  /// @brief To next sensor in the plane
  const Vector3& toNextSensor() const { return m_toNext; }
  /// @brief To next straw layer
  const Vector3& planeNormal() const { return m_planeNorm; }
  /// @brief Measurement's radius
  double driftRadius() const { return m_driftR; }
  /// @brief Measurement's covariance
  const std::array<double, 3>& covariance() const { return m_cov; }
  /// @brief Measurement's drift unceratinty
  double driftUncert() const {
    return std::sqrt(m_cov[toUnderlying(ResidualIdx::bending)]);
  }
  /// @brief Time of record
  double time() const { return m_drifT; }
  /// @brief All measurements are straws
  bool isStraw() const { return m_isStraw; }
  /// @brief Dummy return not used in test
  bool hasTime() const { return false; }
  /// @brief Dummy return not used in test
  bool measuresLoc0() const { return false; }
  /// @brief Dummy return not used in test
  bool measuresLoc1() const { return true; }
  void setRadius(const double r, const double uncertR) {
    m_driftR = abs(r);
    m_cov[toUnderlying(ResidualIdx::bending)] = square(uncertR);
  }
  void setTimeRecord(const double t) { m_drifT = t; }

 private:
  Vector3 m_pos{Vector3::Zero()};
  Vector3 m_wireDir{Vector3::UnitX()};
  Vector3 m_toNext{Vector3::UnitY()};
  Vector3 m_planeNorm{Vector3::UnitZ()};
  double m_driftR{0.};
  std::array<double, 3> m_cov{filledArray<double, 3>(0.)};
  double m_drifT{0.};
  bool m_isStraw{true};
};
static_assert(CompositeSpacePoint<StrawTestPoint>);

class StrawTestCalibrator {
 public:
  /// @brief Choose the coefficient to arrive at a drift time of 750 ns
  ///        for 15 mm
  inline static const double CoeffRtoT = 750._ns / std::pow(15._mm, 1. / 3.);
  inline static const double CoeffTtoR = 1 / std::pow(CoeffRtoT, 3);

  static double calcDriftUncert(const double driftR) {
    return 0.1_mm + 0.15_mm * std::pow(1._mm + std::abs(driftR), -2.);
  }
  static double driftTime(const double r) {
    return CoeffRtoT * std::pow(r, 1. / 3.);
  }
  static double driftRadius(const double t) {
    return CoeffTtoR * std::pow(t, 3);
  }

  static double driftRadius(const CalibrationContext& /*ctx*/,
                            const StrawTestPoint& straw, const double t0) {
    return driftRadius(straw.time() - t0);
  }
  static double driftVelocity(const CalibrationContext& /*ctx*/,
                              const StrawTestPoint& straw, const double t0) {
    return 3 * CoeffTtoR * std::pow(straw.time() - t0, 2);
  }
  static double driftAcceleration(const Acts::CalibrationContext& /*ctx*/,
                                  const StrawTestPoint& straw,
                                  const double t0) {
    return 6 * CoeffTtoR * (straw.time() - t0);
  }
};
static_assert(
    CompositeSpacePointFastCalibrator<StrawTestCalibrator, StrawTestPoint>);

/// @brief Generate a random straight track with a flat distribution in theta & y0
///        Range in theta is [0, 180] degrees in steps of 0.1 excluding the
///        vicinity around 90 degrees.
///          In y0, the generated range is [-500, 500] mm
Line_t generateLine(RandomEngine& engine) {
  using ParIndex = Line_t::ParIndex;
  Line_t::ParamVector linePars{};
  linePars[toUnderlying(ParIndex::x0)] = 0.;
  linePars[toUnderlying(ParIndex::phi)] = 90._degree;
  linePars[toUnderlying(ParIndex::y0)] = uniform_t{-5000., 5000.}(engine);
  linePars[toUnderlying(ParIndex::theta)] =
      uniform_t{0.1_degree, 179.9_degree}(engine);
  if (abs(linePars[toUnderlying(ParIndex::theta)] - 90._degree) < 0.2_degree) {
    return generateLine(engine);
  }
  Line_t line{linePars};
  ACTS_DEBUG("Generated parameters theta: "
             << (linePars[toUnderlying(ParIndex::theta)] / 1._degree)
             << ", y0: " << linePars[toUnderlying(ParIndex::y0)] << " - "
             << toString(line.position()) << " + "
             << toString(line.direction()));
  return line;
}
/// @brief Extrapolate the straight line track through the straw layers to
///        evaluate which tubes were crossed by the track. The straw layers are
///        staggered in the z-direction and each layer expands in y. To
///        estimate, which tubes were crossed, the track is extrapolated to the
///        z-plane below & above the straw wires. Optionally, the true
///        drift-radius can be smeared assuming a Gaussian with a drift-radius
///        dependent uncertainty.
/// @param trajLine: The track to extrapolate
/// @param engine: Random number generator to smear the drift radius
/// @param smearRadius: If true, the drift radius is smeared with a Gaussian
TestStrawCont_t generateStrawCircles(const Line_t& trajLine,
                                     RandomEngine& engine, bool smearRadius) {
  const Vector3 posStaggering{0., std::cos(60._degree), std::sin(60._degree)};
  const Vector3 negStaggering{0., -std::cos(60._degree), std::sin(60._degree)};
  /// Number of tube layers per multilayer
  constexpr std::size_t nLayersPerMl = 8;
  /// Number of overall tubelayers
  constexpr std::size_t nTubeLayers = nLayersPerMl * 2;
  /// Radius of each straw
  constexpr double tubeRadius = 15._mm;
  /// Distance between the first <nLayersPerMl> layers and the second pack
  constexpr double tubeLayerDist = 1.2_m;

  std::array<Vector3, nTubeLayers> tubePositions{
      filledArray<Vector3, nTubeLayers>(Vector3{0., tubeRadius, tubeRadius})};
  /// Fill the positions of the reference tubes 1
  for (std::size_t l = 1; l < nTubeLayers; ++l) {
    const Vector3& layStag{l % 2 == 1 ? posStaggering : negStaggering};
    tubePositions[l] = tubePositions[l - 1] + 2. * tubeRadius * layStag;

    if (l == nLayersPerMl) {
      tubePositions[l] += tubeLayerDist * Vector3::UnitZ();
    }
  }
  /// Print the staggering
  ACTS_DEBUG("##############################################");

  for (std::size_t l = 0; l < nTubeLayers; ++l) {
    ACTS_DEBUG("  *** " << (l + 1) << " - " << toString(tubePositions[l]));
  }
  ACTS_DEBUG("##############################################");

  TestStrawCont_t circles{};
  /// Extrapolate the track to the z-planes of the tubes and determine which
  /// tubes were actually hit
  for (const auto& stag : tubePositions) {
    auto planeExtpLow =
        PlanarHelper::intersectPlane(trajLine.position(), trajLine.direction(),
                                     Vector3::UnitZ(), stag.z() - tubeRadius);
    auto planeExtpHigh =
        PlanarHelper::intersectPlane(trajLine.position(), trajLine.direction(),
                                     Vector3::UnitZ(), stag.z() + tubeRadius);

    ACTS_DEBUG("Extrapolated to plane " << toString(planeExtpLow.position())
                                        << " "
                                        << toString(planeExtpHigh.position()));

    const auto dToFirstLow = static_cast<int>(std::ceil(
        (planeExtpLow.position().y() - stag.y()) / (2. * tubeRadius)));
    const auto dToFirstHigh = static_cast<int>(std::ceil(
        (planeExtpHigh.position().y() - stag.y()) / (2. * tubeRadius)));

    ACTS_DEBUG("Extrapolated to plane "
               << toString(planeExtpLow.position()) << " "
               << toString(planeExtpHigh.position())
               << " Hit tubes: " << dToFirstLow << " " << dToFirstHigh);

    /// Does the track go from left to right or right to left?
    const int dT = dToFirstHigh > dToFirstLow ? 1 : -1;
    /// Loop over the candidate tubes and check each one whether the track
    /// actually crossed them. Then generate the circle and optionally smear the
    /// radius
    for (int tN = dToFirstLow - dT; tN != dToFirstHigh + 2 * dT; tN += dT) {
      const Vector3 tube = stag + 2. * tN * tubeRadius * Vector3::UnitY();
      const double rad = Acts::detail::LineHelper::signedDistance(
          tube, Vector3::UnitX(), trajLine.position(), trajLine.direction());

      if (std::abs(rad) > tubeRadius) {
        continue;
      }
      gauss_t dist{rad, StrawTestCalibrator::calcDriftUncert(rad)};
      const double smearedR = smearRadius ? std::abs(dist(engine)) : rad;
      if (smearedR > tubeRadius) {
        continue;
      }
      circles.emplace_back(std::make_unique<StrawTestPoint>(
          tube, smearedR, StrawTestCalibrator::calcDriftUncert(smearedR)));
      ACTS_DEBUG("Tube position: "
                 << toString(tube) << ", signedRadius: " << rad
                 << ", smearedRadius: " << smearedR << ", uncer: "
                 << StrawTestCalibrator::calcDriftUncert(smearedR));
    }
  }
  ACTS_DEBUG("Track hit in total " << circles.size() << " tubes ");
  return circles;
}

TestStrawCont_t generateStrips(const Line_t& trajLine, RandomEngine& engine) {
  TestStrawCont_t strips{};
  constexpr std::array<double, 16> stripZ{-500._mm, -545_mm, -540_mm, -535._mm,
                                          -250._mm, -245_mm, -240_mm, -235._mm,
                                          335._mm,  340_mm,  345_mm,  350._mm,
                                          435._mm,  440_mm,  445_mm,  400._mm};
  constexpr double stripPitch = 1._cm;
  for (const auto z : stripZ) {
    auto planeExtp = PlanarHelper::intersectPlane(
        trajLine.position(), trajLine.direction(), Vector3::UnitZ(), z);
    Vector3 stripPos = planeExtp.position();
    stripPos[eY] = gauss_t{stripPos[eY], stripPitch}(engine);
    strips.emplace_back(std::make_unique<StrawTestPoint>(stripPos, stripPitch));
  }
  return strips;
}
/// @brief Calculate the overall chi2 of the measurements to the track
/// @param measurements: List of candidate straw measurements
/// @param track: Straight line track
double calcChi2(const TestStrawCont_t& measurements, const Line_t& track) {
  double chi2{0.};
  for (const auto& meas : measurements) {
    const double dist = Acts::detail::LineHelper::signedDistance(
        meas->localPosition(), meas->sensorDirection(), track.position(),
        track.direction());
    ACTS_DEBUG("Distance straw: " << toString(meas->localPosition())
                                  << ",  r: " << meas->driftRadius()
                                  << " - to track: " << abs(dist));

    chi2 += square((abs(dist) - meas->driftRadius()) / meas->driftUncert());
  }
  return chi2;
}

#define DECLARE_BRANCH(dTYPE, NAME) \
  dTYPE NAME{};                     \
  outTree->Branch(#NAME, &NAME);

BOOST_AUTO_TEST_SUITE(SeedingSuite)

void testSimpleStrawFit(RandomEngine& engine, TFile& outFile) {
  auto outTree = std::make_unique<TTree>("StrawFitTree", "FastFitTree");

  DECLARE_BRANCH(double, trueY0);
  DECLARE_BRANCH(double, trueTheta);
  DECLARE_BRANCH(double, recoY0);
  DECLARE_BRANCH(double, recoTheta);
  DECLARE_BRANCH(double, uncertY0);
  DECLARE_BRANCH(double, uncertTheta);
  DECLARE_BRANCH(double, chi2);
  DECLARE_BRANCH(std::size_t, nDoF);
  DECLARE_BRANCH(std::size_t, nIter);

  FastStrawLineFitter::Config cfg{};
  FastStrawLineFitter fastFitter{cfg, getDefaultLogger("FitterNoT0", logLvl)};
  ACTS_INFO("Start simple straw fit test.");
  for (std::size_t n = 0; n < nTrials; ++n) {
    if ((n + 1) % 1000 == 0) {
      ACTS_INFO(" -- processed " << (n + 1) << "/" << nTrials << " events");
    }
    auto track = generateLine(engine);
    auto strawPoints = generateStrawCircles(track, engine, true);
    if (strawPoints.size() < 3) {
      ACTS_WARNING(__func__ << "() - " << __LINE__ << ": -- event: " << n
                            << ", track " << toString(track.position()) << " + "
                            << toString(track.direction())
                            << " did not lead to any valid measurement ");
      continue;
    }
    const std::vector<int> trueDriftSigns =
        CompSpacePointAuxiliaries::strawSigns(track, strawPoints);

    BOOST_CHECK_LE(calcChi2(generateStrawCircles(track, engine, false), track),
                   1.e-12);
    ACTS_DEBUG("True drift signs: " << trueDriftSigns << ", chi2: " << chi2);

    auto fitResult = fastFitter.fit(strawPoints, trueDriftSigns);
    if (!fitResult) {
      continue;
    }
    auto trackPars = track.parameters();

    trueY0 = trackPars[toUnderlying(Line_t::ParIndex::y0)];
    trueTheta = trackPars[toUnderlying(Line_t::ParIndex::theta)];
    /// Calculate the chi2 again
    trackPars[toUnderlying(Line_t::ParIndex::theta)] = (*fitResult).theta;
    trackPars[toUnderlying(Line_t::ParIndex::y0)] = (*fitResult).y0;
    trackPars[toUnderlying(Line_t::ParIndex::phi)] = 90._degree;
    track.updateParameters(trackPars);
    ACTS_DEBUG("Updated parameters: "
               << (trackPars[toUnderlying(Line_t::ParIndex::theta)] / 1._degree)
               << ", y0: " << trackPars[toUnderlying(Line_t::ParIndex::y0)]
               << " -- " << toString(track.position()) << " + "
               << toString(track.direction()));

    const double testChi2 = calcChi2(strawPoints, track);
    ACTS_DEBUG("testChi2: " << testChi2 << ", fit:" << (*fitResult).chi2);

    BOOST_CHECK_LE(abs(testChi2 - (*fitResult).chi2), 1.e-9);
    recoTheta = (*fitResult).theta;
    recoY0 = (*fitResult).y0;
    uncertTheta = (*fitResult).dTheta;
    uncertY0 = (*fitResult).dY0;
    nDoF = (*fitResult).nDoF;
    chi2 = (*fitResult).chi2;
    nIter = (*fitResult).nIter;
    outTree->Fill();
  }
  outFile.WriteObject(outTree.get(), outTree->GetName());
}

void testFitWithT0(RandomEngine& engine, TFile& outFile) {
  auto outTree = std::make_unique<TTree>("StrawFitTreeT0", "FastFitTree");

  DECLARE_BRANCH(double, trueY0);
  DECLARE_BRANCH(double, trueTheta);
  DECLARE_BRANCH(double, trueT0);
  DECLARE_BRANCH(double, recoY0);
  DECLARE_BRANCH(double, recoTheta);
  DECLARE_BRANCH(double, recoT0);
  DECLARE_BRANCH(double, uncertY0);
  DECLARE_BRANCH(double, uncertTheta);
  DECLARE_BRANCH(double, uncertT0);
  DECLARE_BRANCH(double, chi2);
  DECLARE_BRANCH(double, meanSign);
  DECLARE_BRANCH(std::size_t, nDoF);
  DECLARE_BRANCH(std::size_t, nIter);

  FastStrawLineFitter::Config cfg{};
  cfg.maxIter = 500;
  FastStrawLineFitter fastFitter{cfg, getDefaultLogger("FitterWithT0", logLvl)};
  StrawTestCalibrator calibrator{};
  CalibrationContext cctx{};
  ACTS_INFO("Start straw fit with t0 test.");
  for (std::size_t n = 0; n < nTrials; ++n) {
    if ((n + 1) % 1000 == 0) {
      ACTS_INFO(" -- processed " << (n + 1) << "/" << nTrials << " events");
    }
    auto track = generateLine(engine);
    const double timeOffSet = uniform_t{-50._ns, 50._ns}(engine);

    ACTS_DEBUG("Generated time offset: " << inNanoS(timeOffSet) << " [ns]");
    auto strawPoints = generateStrawCircles(track, engine, true);

    if (strawPoints.size() < 4) {
      ACTS_WARNING(__func__ << "() - " << __LINE__ << ": -- event: " << n
                            << ", track " << toString(track.position()) << " + "
                            << toString(track.direction())
                            << " did not lead to any valid measurement ");
      continue;
    }
    BOOST_CHECK_LE(calcChi2(generateStrawCircles(track, engine, false), track),
                   1.e-12);

    /// Fold-in the general offset
    const std::vector<int> trueDriftSigns =
        CompSpacePointAuxiliaries::strawSigns(track, strawPoints);

    ACTS_DEBUG("Straw signs: " << trueDriftSigns);

    for (auto& meas : strawPoints) {
      const double dTime = StrawTestCalibrator::driftTime(meas->driftRadius());
      BOOST_CHECK_CLOSE(StrawTestCalibrator::driftRadius(dTime),
                        meas->driftRadius(), 1.e-12);
      meas->setTimeRecord(dTime + timeOffSet);

      const double updatedR =
          StrawTestCalibrator::driftRadius(dTime + timeOffSet);

      BOOST_CHECK_CLOSE(StrawTestCalibrator::driftRadius(dTime),
                        calibrator.driftRadius(cctx, *meas, timeOffSet), 1.e-3);

      ACTS_DEBUG("Update drift radius of tube "
                 << toString(meas->localPosition()) << " from "
                 << meas->driftRadius() << " to " << updatedR
                 << ", dTime: " << inNanoS(dTime));
      meas->setRadius(updatedR, StrawTestCalibrator::calcDriftUncert(updatedR));

      /// Calculate the numerical derivatives
      constexpr double h = 1.e-8_ns;
      const double numV =
          -(calibrator.driftRadius(cctx, *meas, timeOffSet + h) -
            calibrator.driftRadius(cctx, *meas, timeOffSet - h)) /
          (2. * h);
      BOOST_CHECK_LE(
          abs(numV - calibrator.driftVelocity(cctx, *meas, timeOffSet)) /
              std::max(numV, 1.),
          1.e-3);
    }

    auto result = fastFitter.fit(cctx, calibrator, strawPoints, trueDriftSigns,
                                 std::nullopt);

    if (!result) {
      continue;
    }
    auto linePars = track.parameters();
    /// True parameters
    trueY0 = linePars[toUnderlying(Line_t::ParIndex::y0)];
    trueTheta = linePars[toUnderlying(Line_t::ParIndex::theta)];
    trueT0 = inNanoS(timeOffSet);
    /// Fit parameters
    recoY0 = (*result).y0;
    recoTheta = (*result).theta;
    recoT0 = inNanoS((*result).t0);
    uncertY0 = (*result).dY0;
    uncertTheta = (*result).dTheta;
    uncertT0 = inNanoS((*result).dT0);
    nDoF = (*result).nDoF;
    chi2 = (*result).chi2;
    nIter = (*result).nIter;
    meanSign = {0.};
    std::ranges::for_each(trueDriftSigns,
                          [&meanSign](const int sign) { meanSign += sign; });
    meanSign /= static_cast<double>(trueDriftSigns.size());
    outTree->Fill();
  }
  outFile.WriteObject(outTree.get(), outTree->GetName());
}

void testStripFit(RandomEngine& engine, TFile& outFile) {
  auto outTree = std::make_unique<TTree>("StripFitTree", "FastFitTree");

  DECLARE_BRANCH(double, trueY0);
  DECLARE_BRANCH(double, trueTheta);
  DECLARE_BRANCH(double, recoY0);
  DECLARE_BRANCH(double, recoTheta);
  DECLARE_BRANCH(double, uncertY0);
  DECLARE_BRANCH(double, uncertTheta);
  DECLARE_BRANCH(double, chi2);
  DECLARE_BRANCH(std::size_t, nDoF);
  DECLARE_BRANCH(std::size_t, nIter);

  FastStrawLineFitter::Config cfg{};
  FastStrawLineFitter fastFitter{cfg, getDefaultLogger("StripFitter", logLvl)};
  ACTS_INFO("Start strip fit test.");
  for (std::size_t n = 0; n < nTrials; ++n) {
    if ((n + 1) % 1000 == 0) {
      ACTS_INFO(" -- processed " << (n + 1) << "/" << nTrials << " events");
    }
    auto track = generateLine(engine);
    auto linePars = track.parameters();
    /// True parameters
    trueY0 = linePars[toUnderlying(Line_t::ParIndex::y0)];
    trueTheta = linePars[toUnderlying(Line_t::ParIndex::theta)];

    auto stripPoints = generateStrips(track, engine);

    auto result = fastFitter.fit(stripPoints, ResidualIdx::bending);
    if (!result) {
      continue;
    }
    /// Fit parameters
    recoY0 = (*result).y0;
    recoTheta = (*result).theta;
    uncertY0 = (*result).dY0;
    uncertTheta = (*result).dTheta;
    nDoF = (*result).nDoF;
    nIter = (*result).nIter;
    chi2 = 0.;
    auto pos = recoY0 * Vector3::UnitY();
    auto dir = makeDirectionFromPhiTheta(90._degree, recoTheta);
    for (const auto& strip : stripPoints) {
      chi2 += CompSpacePointAuxiliaries::chi2Term(pos, dir, *strip);
    }

    outTree->Fill();
  }

  outFile.WriteObject(outTree.get(), outTree->GetName());
}

BOOST_AUTO_TEST_CASE(FitterTests) {
  RandomEngine engine{1800};

  std::unique_ptr<TFile> outFile{
      TFile::Open("FastStrawLineFitTest.root", "RECREATE")};

  BOOST_CHECK_EQUAL(outFile->IsZombie(), false);

  testSimpleStrawFit(engine, *outFile);
  testStripFit(engine, *outFile);
  testFitWithT0(engine, *outFile);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests

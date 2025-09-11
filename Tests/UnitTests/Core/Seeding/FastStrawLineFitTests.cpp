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

constexpr std::size_t nTrials = 1;

namespace Acts::Test {

constexpr bool debugMode = true;

ACTS_LOCAL_LOGGER(getDefaultLogger("FastStrawLineFitTests",
                                   Logging::Level::INFO));

class StrawTestPoint;
using TestStrawCont_t = std::vector<std::unique_ptr<StrawTestPoint>>;
using Line_t = CompSpacePointAuxiliaries::Line_t;
using ResidualIdx = FastStrawLineFitter::ResidualIdx;

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

constexpr double inNanoS(const double x) {
  return x / 1._ns;
}

class StrawTestPoint {
 public:
  StrawTestPoint(const Vector3& pos, const double driftR,
                 const double driftRUncert)
      : m_pos{pos}, m_driftR{Acts::abs(driftR)} {
    m_cov[toUnderlying(ResidualIdx::bending)] = Acts::square(driftRUncert);
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
  bool isStraw() const { return true; }
  /// @brief Dummy return not used in test
  bool hasTime() const { return false; }
  /// @brief Dummy return not used in test
  bool measuresLoc0() const { return false; }
  /// @brief Dummy return not used in test
  bool measuresLoc1() const { return false; }
  void setRadius(const double r, const double uncertR) {
    m_driftR = Acts::abs(r);
    m_cov[toUnderlying(ResidualIdx::bending)] = Acts::square(uncertR);
  }
  void setTimeRecord(const double t) { m_drifT = t; }

 private:
  Vector3 m_pos{Vector3::Zero()};
  Vector3 m_wireDir{Vector3::UnitX()};
  Vector3 m_toNext{Vector3::UnitY()};
  Vector3 m_planeNorm{Vector3::UnitZ()};
  double m_driftR{0.};
  std::array<double, 3> m_cov{Acts::filledArray<double, 3>(0.)};
  double m_drifT{0.};
};
static_assert(CompositeSpacePoint<StrawTestPoint>);

class StrawTestCalibrator {
 public:
  /// @brief Choose the coefficient to arrive at a drift time of 750 ns
  ///        for 15 mm
  static constexpr double CoeffRtoT = 750._ns * Acts::pow(15._mm, -2);
  static constexpr double CoeffTtoR = 1. / CoeffRtoT;

  static double calcDriftUncert(const double driftR) {
    return 0.1_mm + 0.15_mm * Acts::pow(1._mm + Acts::abs(driftR), -2);
  }
  static double driftTime(const double r) { return CoeffRtoT * r; }
  static double driftRadius(const double t) {
    return CoeffTtoR * Acts::pow(t, 1);
  }

  static double driftRadius(const Acts::CalibrationContext& /*ctx*/,
                            const StrawTestPoint& straw, const double t0) {
    return driftRadius(straw.time() - t0);
  }
  static double driftVelocity(const Acts::CalibrationContext& /*ctx*/,
                              const StrawTestPoint& straw, const double t0) {
    return CoeffTtoR * Acts::pow(straw.time() - t0, 0);
  }
  static double driftAcceleration(const Acts::CalibrationContext& /*ctx*/,
                                  const StrawTestPoint& /*straw*/,
                                  const double /*t0*/) {
    return 0.;
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
  linePars[toUnderlying(ParIndex::y0)] =
      std::uniform_real_distribution{-5000., 5000.}(engine);
  linePars[toUnderlying(ParIndex::theta)] =
      std::uniform_real_distribution{0.1_degree, 179.9_degree}(engine);
  if (Acts::abs(linePars[toUnderlying(ParIndex::theta)] - 90._degree) <
      0.2_degree) {
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
    auto planeExtpLow = Acts::PlanarHelper::intersectPlane(
        trajLine.position(), trajLine.direction(), Vector3::UnitZ(),
        stag.z() - tubeRadius);
    auto planeExtpHigh = Acts::PlanarHelper::intersectPlane(
        trajLine.position(), trajLine.direction(), Vector3::UnitZ(),
        stag.z() + tubeRadius);

    ACTS_DEBUG("Extrapolated to plane " << toString(planeExtpLow.position())
                                        << " "
                                        << toString(planeExtpHigh.position()));

    const auto dToFirstLow = static_cast<int>(std::ceil(
        (planeExtpLow.position().y() - stag.y()) / (2. * tubeRadius)));
    const auto dToFirstHigh = static_cast<int>(std::ceil(
        (planeExtpHigh.position().y() - stag.y()) / (2. * tubeRadius)));
    /// Does the track go from left to right or right to left?
    const int dT = dToFirstHigh > dToFirstLow ? 1 : -1;
    /// Loop over the candidate tubes and check each one whether the track
    /// actually crossed them. Then generate the circle and optionally smear the
    /// radius
    for (int tN = dToFirstLow - dT; tN != dToFirstHigh + 2 * dT; tN += dT) {
      const Vector3 tube = stag + 2. * tN * tubeRadius * Vector3::UnitY();
      const double rad = Acts::detail::LineHelper::signedDistance(
          tube, Vector3::UnitX(), trajLine.position(), trajLine.direction());
      ACTS_DEBUG("Tube position: " << toString(tube) << ", radius: " << rad);

      if (std::abs(rad) > tubeRadius) {
        continue;
      }
      std::normal_distribution<> dist{
          rad, StrawTestCalibrator::calcDriftUncert(rad)};
      const double smearedR = smearRadius ? std::abs(dist(engine)) : rad;
      if (smearedR > tubeRadius) {
        continue;
      }
      circles.emplace_back(std::make_unique<StrawTestPoint>(
          tube, smearedR, StrawTestCalibrator::calcDriftUncert(smearedR)));
    }
  }
  ACTS_DEBUG("Track hit in total " << circles.size() << " tubes ");
  return circles;
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
                                  << " - to track: " << Acts::abs(dist));

    chi2 += Acts::square((Acts::abs(dist) - meas->driftRadius()) /
                         meas->driftUncert());
  }
  return chi2;
}

BOOST_AUTO_TEST_SUITE(FastStrawLineFitTests)

BOOST_AUTO_TEST_CASE(SimpleLineFit) {
  RandomEngine engine{1419};

  std::unique_ptr<TFile> outFile{};
  std::unique_ptr<TTree> outTree{};
  double trueY0{0.};
  double trueTheta{0.};
  double fitY0{0.};
  double fitTheta{0.};
  double fitdY0{0.};
  double fitdTheta{0.};
  double chi2{0.};
  std::size_t nDoF{0u};
  std::size_t nIter{0u};
  if (debugMode) {
    outFile.reset(TFile::Open("FastStrawLineFitTest.root", "RECREATE"));
    BOOST_CHECK_EQUAL(outFile->IsZombie(), false);
    outTree = std::make_unique<TTree>("FastFitTree", "FastFitTree");
    outTree->Branch("trueY0", &trueY0);
    outTree->Branch("trueTheta", &trueTheta);
    outTree->Branch("fitY0", &fitY0);
    outTree->Branch("fitTheta", &fitTheta);
    outTree->Branch("errY0", &fitdY0);
    outTree->Branch("errTheta", &fitdTheta);
    outTree->Branch("chi2", &chi2);
    outTree->Branch("nDoF", &nDoF);
    outTree->Branch("nIter", &nIter);
  }

  FastStrawLineFitter::Config cfg{};
  FastStrawLineFitter fastFitter{cfg};
  for (std::size_t n = 0; n < nTrials; ++n) {
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

    BOOST_CHECK_LE(Acts::abs(testChi2 - (*fitResult).chi2), 1.e-9);
    if (debugMode) {
      fitTheta = (*fitResult).theta;
      fitY0 = (*fitResult).y0;
      fitdTheta = (*fitResult).dTheta;
      fitdY0 = (*fitResult).dY0;
      nDoF = (*fitResult).nDoF;
      chi2 = (*fitResult).chi2;
      nIter = (*fitResult).nIter;
      outTree->Fill();
    }
  }
  if (debugMode) {
    outFile->WriteObject(outTree.get(), outTree->GetName());
    outTree.reset();
  }
}

BOOST_AUTO_TEST_CASE(LineFitWithT0) {
  RandomEngine engine{47110};

  std::unique_ptr<TFile> outFile{};
  std::unique_ptr<TTree> outTree{};
  double trueY0{0.};
  double trueTheta{0.};
  double trueT0{0.};
  double fitY0{0.};
  double fitTheta{0.};
  double fitT0{0.};
  double fitdY0{0.};
  double fitdTheta{0.};
  double fitDT0{0.};
  double chi2{0.};
  double meanSign{0.};
  std::size_t nDoF{0u};
  std::size_t nIter{0u};
  if (debugMode) {
    outFile.reset(TFile::Open("FastStrawLineFitTestT0.root", "RECREATE"));
    BOOST_CHECK_EQUAL(outFile->IsZombie(), false);
    outTree = std::make_unique<TTree>("FastFitTree", "FastFitTree");
    outTree->Branch("trueY0", &trueY0);
    outTree->Branch("trueTheta", &trueTheta);
    outTree->Branch("trueT0", &trueT0);
    outTree->Branch("fitY0", &fitY0);
    outTree->Branch("fitTheta", &fitTheta);
    outTree->Branch("fitT0", &fitT0);

    outTree->Branch("errY0", &fitdY0);
    outTree->Branch("errTheta", &fitdTheta);
    outTree->Branch("errT0", &fitDT0);

    outTree->Branch("chi2", &chi2);
    outTree->Branch("nDoF", &nDoF);
    outTree->Branch("nIter", &nIter);
  }

  FastStrawLineFitter::Config cfg{};
  cfg.maxIter = 500;
  FastStrawLineFitter fastFitter{cfg};
  StrawTestCalibrator calibrator{};
  Acts::CalibrationContext cctx{};
  for (std::size_t n = 0; n < nTrials; ++n) {
    auto track = generateLine(engine);
    const double timeOffSet =
        std::uniform_real_distribution{-50._ns, 50._ns}(engine);

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

      const double updatedR =
          StrawTestCalibrator::driftRadius(dTime + timeOffSet);

      meas->setTimeRecord(dTime + timeOffSet);
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
          Acts::abs(numV - calibrator.driftVelocity(cctx, *meas, timeOffSet)) /
              std::max(numV, 1.),
          1.e-3);
    }

    auto result = fastFitter.fit(cctx, calibrator, strawPoints, trueDriftSigns,
                                 std::nullopt);

    if (!result) {
      continue;
    }
    auto linePars = track.parameters();
    if (debugMode) {
      /// True parameters
      trueY0 = linePars[toUnderlying(Line_t::ParIndex::y0)];
      trueTheta = linePars[toUnderlying(Line_t::ParIndex::theta)];
      trueT0 = inNanoS(timeOffSet);
      /// Fit parameters
      fitY0 = (*result).y0;
      fitTheta = (*result).theta;
      fitT0 = inNanoS((*result).t0);
      fitdY0 = (*result).dY0;
      fitdTheta = (*result).dTheta;
      fitDT0 = inNanoS((*result).dT0);
      nDoF = (*result).nDoF;
      chi2 = (*result).chi2;
      nIter = (*result).nIter;
      meanSign = {0.};
      std::ranges::for_each(trueDriftSigns,
                            [&meanSign](const int sign) { meanSign += sign; });
      meanSign /= static_cast<double>(trueDriftSigns.size());
      outTree->Fill();
    }
  }
  if (debugMode) {
    outFile->WriteObject(outTree.get(), outTree->GetName());
    outTree.reset();
  }
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Acts::Test

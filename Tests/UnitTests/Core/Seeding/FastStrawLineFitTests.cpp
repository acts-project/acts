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

constexpr std::uint32_t nTrials = 10000;

namespace Acts::Test {

constexpr bool debugMode = true;
constexpr bool print = false;
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

class StrawTestPoint {
 public:
  StrawTestPoint(const Vector3& pos, const double driftR,
                 const double driftRUncert)
      : m_pos{pos}, m_driftR{Acts::abs(driftR)} {
    m_cov[toUnderlying(ResidualIdx::bending)] = Acts::pow(driftRUncert, 2);
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
  /// @brief Dummy return not used in test
  double time() const { return 0.; }
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
    m_cov[toUnderlying(ResidualIdx::bending)] = Acts::pow(uncertR, 2);
  }

 private:
  Vector3 m_pos{Vector3::Zero()};
  Vector3 m_wireDir{Vector3::UnitX()};
  Vector3 m_toNext{Vector3::UnitY()};
  Vector3 m_planeNorm{Vector3::UnitZ()};
  double m_driftR{0.};
  std::array<double, 3> m_cov{Acts::filledArray<double, 3>(0.)};
};
static_assert(CompositeSpacePoint<StrawTestPoint>);

class StrawTestCalibrator {
 public:
  static double driftVelocity(const Acts::CalibrationContext& /*ctx*/,
                              const StrawTestPoint& straw) {
    constexpr double A = 1. / (750._ns * Acts::pow(15._mm, -2));
    const double t = driftTime(straw.driftRadius());
    return A / (2. * std::sqrt(t * A));
  }
  static double driftAcceleration(const Acts::CalibrationContext& /*ctx*/,
                                  const StrawTestPoint& straw) {
    constexpr double A = 1. / (750._ns * Acts::pow(15._mm, -2));
    const double t = driftTime(straw.driftRadius());

    return A / (2. * driftRadius(t * A));
  }

  static constexpr double calcDriftUncert(const double driftR) {
    return 0.1_mm + 0.15_mm * Acts::pow(1._mm + Acts::abs(driftR), -2);
  }
  static constexpr double driftTime(const double r) {
    constexpr double A = 750._ns * Acts::pow(15._mm, -2);
    return A * Acts::pow(r, 2);
  }
  static double driftRadius(const double t) {
    constexpr double A = 1. / (750._ns * Acts::pow(15._mm, -2));
    return std::sqrt(t * A);
  }
};
static_assert(
    CompositeSpacePointCalibrator<StrawTestCalibrator, StrawTestPoint>);

Line_t generateLine(RandomEngine& engine) {
  using ParIndex = Line_t::ParIndex;
  Line_t::ParamVector linePars{};
  linePars[toUnderlying(ParIndex::x0)] = 0.;
  linePars[toUnderlying(ParIndex::phi)] = 90._degree;
  linePars[toUnderlying(ParIndex::y0)] = (engine() % 10000 - 5000.) / 10.;
  constexpr unsigned maxAngle = 180;
  linePars[toUnderlying(ParIndex::theta)] =
      (engine() % (10 * maxAngle)) * 0.1_degree;
  Line_t line{};
  line.updateParameters(linePars);
  if (Acts::abs(linePars[toUnderlying(ParIndex::theta)] - 90._degree) <
      0.2_degree) {
    return generateLine(engine);
  }
  if constexpr (print) {
    std::cout << "Generated parameters theta: "
              << (linePars[toUnderlying(ParIndex::theta)] / 1._degree)
              << ", y0: " << linePars[toUnderlying(ParIndex::y0)] << " - "
              << toString(line.position()) << " + "
              << toString(line.direction()) << std::endl;
  }
  return line;
}

TestStrawCont_t generateStrawCircles(const Line_t& trajLine,
                                     RandomEngine& engine, bool smearRadius) {
  const Vector3 posStaggering{0., std::cos(60._degree), std::sin(60._degree)};
  const Vector3 negStaggering{0., -std::cos(60._degree), std::sin(60._degree)};
  /// Number of tube layers per multilayer
  constexpr std::uint32_t nLayersPerMl = 8;
  /// Number of overall tubelayers
  constexpr std::uint32_t nTubeLayers = nLayersPerMl * 2;
  constexpr double tubeRadius = 15._mm;
  constexpr double tubeLayerDist = 1.2_m;

  std::array<Vector3, nTubeLayers> tubePositions{
      filledArray<Vector3, nTubeLayers>(Vector3{0., tubeRadius, tubeRadius})};

  for (std::uint32_t l = 1; l < nTubeLayers; ++l) {
    const Vector3& layStag{l % 2 == 1 ? posStaggering : negStaggering};
    tubePositions[l] = tubePositions[l - 1] + 2. * tubeRadius * layStag;

    if (l == nLayersPerMl) {
      tubePositions[l] += tubeLayerDist * Vector3::UnitZ();
    }
  }
  /// Print the staggering
  if constexpr (print) {
    std::cout << "##############################################" << std::endl;
    for (std::uint32_t l = 0; l < nTubeLayers; ++l) {
      std::cout << " *** " << (l + 1) << " - " << toString(tubePositions[l])
                << std::endl;
    }
    std::cout << "##############################################" << std::endl;
  }
  TestStrawCont_t circles{};
  for (const auto& stag : tubePositions) {
    auto planeExtpLow = Acts::PlanarHelper::intersectPlane(
        trajLine.position(), trajLine.direction(), Vector3::UnitZ(),
        stag.z() - tubeRadius);
    auto planeExtpHigh = Acts::PlanarHelper::intersectPlane(
        trajLine.position(), trajLine.direction(), Vector3::UnitZ(),
        stag.z() + tubeRadius);
    if constexpr (print) {
      std::cout << "extrapolated to plane " << toString(planeExtpLow.position())
                << " " << toString(planeExtpHigh.position()) << std::endl;
    }
    const auto dToFirstLow = static_cast<std::int32_t>(std::ceil(
        (planeExtpLow.position().y() - stag.y()) / (2. * tubeRadius)));
    const auto dToFirstHigh = static_cast<std::int32_t>(std::ceil(
        (planeExtpHigh.position().y() - stag.y()) / (2. * tubeRadius)));

    const std::int32_t dT = dToFirstHigh > dToFirstLow ? 1 : -1;
    for (std::int32_t tN = dToFirstLow - dT; tN != dToFirstHigh + 2 * dT;
         tN += dT) {
      const Vector3 tube = stag + 2. * tN * tubeRadius * Vector3::UnitY();
      const double rad = Acts::detail::LineHelper::signedDistance(
          tube, Vector3::UnitX(), trajLine.position(), trajLine.direction());
      if constexpr (print) {
        std::cout << "Tube position: " << toString(tube) << ", radius: " << rad
                  << std::endl;
      }
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
  if constexpr (print) {
    std::cout << "Track hit in total " << circles.size() << " tubes "
              << std::endl;
  }
  return circles;
}

double calcChi2(const TestStrawCont_t& measurements, const Line_t& track) {
  double chi2{0.};
  for (const auto& meas : measurements) {
    const double dist = Acts::detail::LineHelper::signedDistance(
        meas->localPosition(), meas->sensorDirection(), track.position(),
        track.direction());
    if constexpr (print) {
      std::cout << "calcChi2() - Distance straw: "
                << toString(meas->localPosition())
                << ",  r: " << meas->driftRadius()
                << " - to track: " << Acts::abs(dist) << std::endl;
    }
    chi2 += Acts::pow(
        (Acts::abs(dist) - meas->driftRadius()) / meas->driftUncert(), 2);
  }
  return chi2;
}

BOOST_AUTO_TEST_SUITE(FastStrawLineFitTests)

BOOST_AUTO_TEST_CASE(SimpleLineFit) {
  RandomEngine engine{1419};
  return;

  std::unique_ptr<TFile> outFile{};
  std::unique_ptr<TTree> outTree{};
  double trueY0{0.}, trueTheta{0.}, fitY0{0.}, fitTheta{0.};
  double fitdY0{0.}, fitdTheta{0.}, chi2{0.};
  std::uint32_t nDoF{0u}, nIter{0u};
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
  for (std::uint32_t n = 0; n < nTrials; ++n) {
    auto track = generateLine(engine);
    auto strawPoints = generateStrawCircles(track, engine, true);
    if (strawPoints.size() < 3) {
      std::cout << "WARNING -- event: " << n << ", track "
                << toString(track.position()) << " + "
                << toString(track.direction())
                << " did not lead to any valid measurement " << std::endl;
      continue;
    }
    std::vector<std::int32_t> trueDriftSigns{};
    trueDriftSigns.reserve(strawPoints.size());
    for (const auto& meas : strawPoints) {
      trueDriftSigns.push_back(
          CompSpacePointAuxiliaries::strawSign(track, *meas));
    }
    BOOST_CHECK_LE(calcChi2(generateStrawCircles(track, engine, false), track),
                   1.e-12);
    if constexpr (print) {
      std::cout << "True drift signs: " << trueDriftSigns << ", chi2: " << chi2
                << std::endl;
    }
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
    if constexpr (print) {
      std::cout << "Updated parameters: "
                << (trackPars[toUnderlying(Line_t::ParIndex::theta)] /
                    1._degree)
                << ", y0: " << trackPars[toUnderlying(Line_t::ParIndex::y0)]
                << " -- " << toString(track.position()) << " + "
                << toString(track.direction()) << std::endl;
    }

    const double testChi2 = calcChi2(strawPoints, track);
    if constexpr (print) {
      std::cout << "testChi2: " << testChi2 << ", fit:" << (*fitResult).chi2
                << std::endl;
    }
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

  FastStrawLineFitter::Config cfg{};
  FastStrawLineFitter fastFitter{cfg};
  StrawTestCalibrator calibrator{};
  Acts::CalibrationContext ctx{};
  for (std::uint32_t n = 0; n < nTrials; ++n) {
    auto track = generateLine(engine);
    const double timeOffSet = 5._ns + (engine() % 100) * 1._ns;

    auto strawPoints = generateStrawCircles(track, engine, true);
    if (strawPoints.size() < 4) {
      std::cout << "WARNING -- event: " << n << ", track "
                << toString(track.position()) << " + "
                << toString(track.direction())
                << " did not lead to any valid measurement " << std::endl;
      continue;
    }
    /// Fold-in the general offset
    std::vector<std::int32_t> trueDriftSigns{};
    trueDriftSigns.reserve(strawPoints.size());

    for (auto& meas : strawPoints) {
      const double dTime = StrawTestCalibrator::driftTime(meas->driftRadius());
      BOOST_CHECK_CLOSE(StrawTestCalibrator::driftRadius(dTime),
                        meas->driftRadius(), 1.e-12);
      const double updatedR =
          StrawTestCalibrator::driftRadius(dTime + timeOffSet);
      meas->setRadius(updatedR, StrawTestCalibrator::calcDriftUncert(updatedR));
      trueDriftSigns.push_back(
          CompSpacePointAuxiliaries::strawSign(track, *meas));
    }
    auto result = fastFitter.fit(ctx, calibrator, strawPoints, trueDriftSigns);
    /// Bail out
    break;

    if (!result) {
      continue;
    }
   
  }
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Acts::Test

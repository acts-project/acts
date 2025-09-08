// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/CompositeSpacePointLineFitter.hpp"

#include <random>

using namespace Acts;
using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;

using RandomEngine = std::mt19937;
using uniform = std::uniform_real_distribution<double>;
using Line_t = CompSpacePointAuxiliaries::Line_t;
constexpr auto logLvl = Acts::Logging::Level::INFO;

namespace Acts::Test {

class FitTestSpacePoint {
 public:
  /// @brief Constructor for straw wires
  /// @param pos: Position of the wire
  /// @param driftR: Straw drift radius
  /// @param driftRUncert: Uncertainty on the drift radius uncertainty
  /// @param mLoc0: Flag toggling whether the position along the wire is also provided
  FitTestSpacePoint(const Vector3& pos, const double driftR,
                    const double driftRUncert, const bool mLoc0 = false)
      : m_position{pos},
        m_driftR{driftR},
        m_measLoc0{mLoc0} {m_covariance[toUnderlying()] = }

        /// @brief Constructor for strip measurements
        FitTestSpacePoint(const Vector3& stripPos, const Vector3& stripDir,
                          const Vector3& toNext, const double uncertLoc0,
                          const double uncertLoc1)
      : m_position{pos},
        m_sensorDir{stripDir},
        m_toNextSen{toNext},
        m_measLoc0{uncertLoc0 > 0},
        m_measLoc1{uncertLoc1 > 0} {}

  const Vector3& localPosition() const { return m_position; }
  /// @brief Wire direction
  const Vector3& sensorDirection() const { return m_sensorDir; }
  /// @brief To next sensor in the plane
  const Vector3& toNextSensor() const { return m_toNextSen; }
  /// @brief To next straw layer
  const Vector3& planeNormal() const { return m_planeNorm; }
  /// @brief Measurement's radius
  double driftRadius() const { return m_driftR.value_or(0.); }
  /// @brief Measurement's covariance
  const std::array<double, 3>& covariance() const { return m_covariance; }
  /// @brief Time of record
  double time() const { return m_time.value_or(0.); }
  /// @brief All measurements are straws
  bool isStraw() const { return m_driftR.has_value(); }
  /// @brief Check whether the space point has a time value
  bool hasTime() const { return m_time.has_value(); }
  /// @brief Check whether the space point measures the non-bending direction
  bool measuresLoc0() const { return m_measLoc0; }
  /// @brief Check whether the space point measures the bending direction
  bool measuresLoc1() const { return m_measLoc1 || isStraw(); }
  /// @brief Sets the straw tube's drift radius
  void updateDriftR(const double updatedR) { m_driftR = updatedR; }

 private:
  Vector3 m_position{Vector3::Zero()};
  Vector3 m_sensorDir{Vector3::UnitX()};
  Vector3 m_toNextSen{Vector3::UnitY()};
  Vector3 m_planeNorm{m_sensorDir.cross(m_toNextSen).normalized()};
  std::optional<double> m_driftR{std::nullopt};
  std::optional<double> m_time{std::nullopt};
  std::array<double, 3> m_covariance{filledArray<double, 3>(0)};
  bool m_measLoc0{false};
  bool m_measLoc1{false};
};

static_assert(CompositeSpacePoint<FitTestSpacePoint>);

using Container_t = std::vector<std::shared_ptr<FitTestSpacePoint>>;

static_assert(CompositeSpacePointContainer<Container_t>);
class SpCalibrator {
 public:
  Container_t calibrate(const Acts::CalibrationContext& /*ctx*/,
                        const Vector3& trackPos, const Vector3& trackDir,
                        const double timeOffSet,
                        const Container_t& uncalibCont) const {
    return uncalibCont;
  }

  void updateSigns(const Vector3& trackPos, const Vector3& trackDir,
                   Container_t& measurements) const {
    auto signs =
        CompSpacePointAuxiliaries::strawSigns(trackPos, trackDir, measurements);
    for (std::size_t s = 0; s < signs.size(); ++s) {
      if (measurements[s]->isStraw()) {
        measurements[s]->updateDriftR(
            Acts::abs(measurements[s]->driftRadius() * signs[s]));
      }
    }
  }
};

static_assert(
    CompositeSpacePointCalibrator<SpCalibrator, Container_t, Container_t>);

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
Container_t generateMeasurements(const Line_t& trajLine, const double t0) {
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

ACTS_LOCAL_LOGGER(getDefaultLogger("StrawLineFitterTest", logLvl));

BOOST_AUTO_TEST_SUITE(FastStrawLineFitTests)

BOOST_AUTO_TEST_CASE(SimpleLineFit) {
  RandomEngine engine{1503};

  using FitOpts_t =
      CompositeSpacePointLineFitter::FitOptions<Container_t, SpCalibrator>;

  CompositeSpacePointLineFitter::Config cfg{};

  CompositeSpacePointLineFitter fitter{cfg};

  FitOpts_t fitOpts{};
  auto result = fitter.fit(std::move(fitOpts));
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Acts::Test

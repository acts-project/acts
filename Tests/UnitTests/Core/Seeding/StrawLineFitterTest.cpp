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
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <random>

using namespace Acts;
using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;
using namespace Acts::UnitLiterals;
using namespace Acts::detail::LineHelper;
using namespace Acts::PlanarHelper;
using namespace Acts::VectorHelpers;

using RandomEngine = std::mt19937;
using uniform = std::uniform_real_distribution<double>;
using normal_t = std::normal_distribution<double>;
using Line_t = CompSpacePointAuxiliaries::Line_t;
using ResidualIdx = CompSpacePointAuxiliaries::ResidualIdx;
constexpr auto logLvl = Acts::Logging::Level::INFO;

ACTS_LOCAL_LOGGER(getDefaultLogger("StrawLineFitterTest", logLvl));

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
      : m_position{pos}, m_driftR{driftR}, m_measLoc0{mLoc0} {
    m_covariance[toUnderlying(ResidualIdx::bending)] =
        Acts::square(driftRUncert);
  }

  /// @brief Constructor for strip measurements
  FitTestSpacePoint(const Vector3& stripPos, const Vector3& stripDir,
                    const Vector3& toNext, const double uncertLoc0,
                    const double uncertLoc1)
      : m_position{stripPos},
        m_sensorDir{stripDir},
        m_toNextSen{toNext},
        m_measLoc0{uncertLoc0 > 0.},
        m_measLoc1{uncertLoc1 > 0.} {
    using enum ResidualIdx;
    m_covariance[toUnderlying(nonBending)] = Acts::square(uncertLoc0);
    m_covariance[toUnderlying(bending)] = Acts::square(uncertLoc1);
  }
  /// @brief Position of the space point
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
  /// @brief Updates the position of the space point
  void updatePosition(const Vector3& newPos) { m_position = newPos; }

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
  static double driftUncert(const double r) {
    return 0.2_mm / (1._mm + Acts::square(r)) + 0.1_mm;
  }
  /// @brief Calibrate a set of straw measurements using the best known estimate on a straight line track
  /// @param ctx: Calibration context (Needed by conept interface)
  /// @param trackPos: Position of the track at z=0.
  /// @param trackDir: Direction of the track in the local frame
  /// @param timeOffSet: Offset in the time of arrival (To be implemented)
  /// @param uncalibCont: Uncalibrated composite space point container
  Container_t calibrate(const Acts::CalibrationContext& /*ctx*/,
                        const Vector3& trackPos, const Vector3& trackDir,
                        const double /*timeOffSet*/,
                        const Container_t& uncalibCont) const {
    Container_t calibMeas{};
    for (const auto& sp : uncalibCont) {
      if (!sp->measuresLoc0() || !sp->measuresLoc1()) {
        /// Estimate the best position along the sensor
        auto bestPos = lineIntersect(trackPos, trackDir, sp->localPosition(),
                                     sp->sensorDirection());
        sp->updatePosition(bestPos.position());
      }
      if (sp->isStraw()) {
        sp->updateDriftR(Acts::abs(sp->driftRadius()));
      }
    }
    return uncalibCont;
  }
  /// @brief Updates the sign of the Straw's drift radii indicating that they are on the left (-1)
  ///        or right side (+1) of the track line
  void updateSigns(const Vector3& trackPos, const Vector3& trackDir,
                   Container_t& measurements) const {
    auto signs =
        CompSpacePointAuxiliaries::strawSigns(trackPos, trackDir, measurements);
    /// The signs have the same size as the measurement container
    for (std::size_t s = 0; s < signs.size(); ++s) {
      /// Take care to not turn strips into straws by accident
      if (measurements[s]->isStraw()) {
        measurements[s]->updateDriftR(
            Acts::abs(measurements[s]->driftRadius()) * signs[s]);
      }
    }
  }
};
/// Ensure that the Test space point calibrator satisfies the calibrator concept
static_assert(
    CompositeSpacePointCalibrator<SpCalibrator, Container_t, Container_t>);

/// @brief Generates a random straight line
/// @param engine Random number sequence to draw the parameters from
/// @return A Line object instantiated with the generated parameters
Line_t generateLine(RandomEngine& engine) {
  using ParIndex = Line_t::ParIndex;
  Line_t::ParamVector linePars{};
  linePars[toUnderlying(ParIndex::x0)] =
      std::uniform_real_distribution{-120_degree, 120_degree}(engine);
  linePars[toUnderlying(ParIndex::phi)] =
      std::uniform_real_distribution{-5000., 5000.}(engine);
  linePars[toUnderlying(ParIndex::y0)] =
      std::uniform_real_distribution{-5000., 5000.}(engine);
  linePars[toUnderlying(ParIndex::theta)] =
      std::uniform_real_distribution{0.1_degree, 179.9_degree}(engine);
  if (Acts::abs(linePars[toUnderlying(ParIndex::theta)] - 90._degree) <
      0.2_degree) {
    return generateLine(engine);
  }
  Line_t line{linePars};
  ACTS_DEBUG(
      "Generated parameters -"
      << std::format("theta: {:.2f}, phi: {:.2f}, y0: {:.1f}, x0: {:.1f}",
                     (linePars[toUnderlying(ParIndex::theta)] / 1._degree),
                     (linePars[toUnderlying(ParIndex::phi)] / 1._degree),
                     linePars[toUnderlying(ParIndex::y0)],
                     linePars[toUnderlying(ParIndex::x0)])
      << " --> " << toString(line.position()) << " + "
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
Container_t generateMeasurements(const Line_t& trajLine, const double /*t0*/,
                                 RandomEngine& engine,
                                 const bool smearRadius = true) {
  /// Direction vector to go a positive step in the tube honeycomb grid
  const Vector3 posStaggering{0., std::cos(60._degree), std::sin(60._degree)};
  /// Direction vector to go a negative step in the tube honeycomb grid
  const Vector3 negStaggering{0., -std::cos(60._degree), std::sin(60._degree)};
  /// Number of tube layers per multilayer
  constexpr std::size_t nLayersPerMl = 8;
  /// Number of overall tubelayers
  constexpr std::size_t nTubeLayers = nLayersPerMl * 2;
  /// Radius of each straw
  constexpr double tubeRadius = 15._mm;
  /// Distance between the first <nLayersPerMl> layers and the second pack
  constexpr double tubeLayerDist = 1.2_m;
  /// Distance between the tube multilayers and to the first strip layer
  constexpr double tubeStripDist = 30._cm;
  /// Distance between two strip layers
  constexpr double stripLayDist = 0.5_cm;
  /// Strip pitch in x0 direction
  constexpr double stripPitchLoc0 = 4._cm;
  /// Strip pitch in y0 direction
  constexpr double stripPitchLoc1 = 2._cm;
  /// Number of strip layers on each side
  constexpr std::size_t nStripLay = 4;

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

  Container_t measurements{};
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
      const double smearedR =
          smearRadius
              ? std::abs(normal_t{rad, SpCalibrator::driftUncert(rad)}(engine))
              : rad;
      if (smearedR > tubeRadius) {
        continue;
      }
      measurements.emplace_back(std::make_unique<FitTestSpacePoint>(
          tube, smearedR, SpCalibrator::driftUncert(smearedR)));
    }
  }
  for (std::size_t sL = 0; sL < nStripLay; ++sL) {
    const double distInZ = tubeStripDist + sL * stripLayDist;
    const double planeLow = tubePositions.front().z() - distInZ;
    const double planeHigh = tubePositions.back().z() + distInZ;
    for (const double plane : {planeLow, planeHigh}) {
      const auto extp = intersectPlane(
          trajLine.position(), trajLine.direction(), Vector3::UnitZ(), plane);
      Vector3 extpPos = extp.position();
      extpPos[0] = normal_t{extpPos[0], stripPitchLoc0}(engine);
      extpPos[1] = 0.;
      measurements.emplace_back(std::make_unique<FitTestSpacePoint>(
          extpPos, Vector3::UnitY(), Vector3::UnitX(), stripPitchLoc0, 0.));
      extpPos = extp.position();
      extpPos[1] = normal_t{extpPos[1], stripPitchLoc1}(engine);
      extpPos[0] = 0.;
      measurements.emplace_back(std::make_unique<FitTestSpacePoint>(
          extpPos, Vector3::UnitX(), Vector3::UnitY(), 0., stripPitchLoc1));
    }
  }
  ACTS_DEBUG("Track hit in total " << measurements.size() << " sensors.");
  return measurements;
}

BOOST_AUTO_TEST_SUITE(FastStrawLineFitTests)

BOOST_AUTO_TEST_CASE(SimpleLineFit) {
  RandomEngine engine{1503};

  using FitOpts_t =
      CompositeSpacePointLineFitter::FitOptions<Container_t, SpCalibrator>;

  CompositeSpacePointLineFitter::Config cfg{};

  CompositeSpacePointLineFitter fitter{cfg,
                                       getDefaultLogger("LineFitter", logLvl)};

  auto calibrator = std::make_unique<SpCalibrator>();
  for (std::size_t evt = 0; evt < 10; ++evt) {
    const auto line = generateLine(engine);
    const double t0 = uniform{-50._ns, 50._ns}(engine);

    FitOpts_t fitOpts{};
    fitOpts.calibrator = calibrator.get();
    fitOpts.measurements = generateMeasurements(line, t0, engine);

    double tanPhi{0.};
    double tanTheta{0.};
    /// Setup the seed parameters in x0 && phi
    {
      auto firstPhi = std::ranges::find_if(
          fitOpts.measurements,
          [](const auto& sp) { return !sp->isStraw() && sp->measuresLoc0(); });
      auto lastPhi = std::ranges::find_if(
          std::ranges::reverse_view(fitOpts.measurements),
          [](const auto& sp) { return !sp->isStraw() && sp->measuresLoc0(); });
      const Vector firstToLastPhi =
          (**lastPhi).localPosition() - (**firstPhi).localPosition();
      tanPhi = firstToLastPhi.x() / firstToLastPhi.z();
      /// -> x = tanPhi * z + x_{0} ->
      fitOpts.startParameters[toUnderlying(FitParIndex::x0)] =
          (**lastPhi).localPosition().x() -
          (**lastPhi).localPosition().z() * tanPhi;
    }
    /// Setup the seed parameters in y0 && theta
    {
      auto firstTube = std::ranges::find_if(
          fitOpts.measurements, [](const auto& sp) { return sp->isStraw(); });
      auto lastTube =
          std::ranges::find_if(std::ranges::reverse_view(fitOpts.measurements),
                               [](const auto& sp) { return sp->isStraw(); });
      const int signFirst =
          CompSpacePointAuxiliaries::strawSign(line, **firstTube);
      const int signLast =
          CompSpacePointAuxiliaries::strawSign(line, **lastTube);
    }

    const Vector3 seedDir = makeDirectionFromAxisTangents(tanPhi, tanTheta);
    fitOpts.startParameters[toUnderlying(FitParIndex::theta)] = theta(seedDir);
    fitOpts.startParameters[toUnderlying(FitParIndex::phi)] = phi(seedDir);

    //
    auto result = fitter.fit(std::move(fitOpts));
  }
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Acts::Test

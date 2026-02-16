// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/CompositeSpacePointLineFitter.hpp"
#include "Acts/Seeding/CompositeSpacePointLineSeeder.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/detail/Polynomials.hpp"

#include <random>
#include <ranges>

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace Acts::VectorHelpers;
using namespace Acts::Experimental::detail;
using namespace Acts::Experimental;
using namespace Acts::detail::LineHelper;
using namespace Acts::PlanarHelper;

using RandomEngine = std::mt19937;

using ResidualIdx = CompSpacePointAuxiliaries::ResidualIdx;
using Line_t = CompSpacePointAuxiliaries::Line_t;
using normal_t = std::normal_distribution<double>;
using uniform = std::uniform_real_distribution<double>;
using FitParIndex = CompSpacePointAuxiliaries::FitParIndex;
using ParamVec_t = CompositeSpacePointLineFitter::ParamVec_t;

namespace ActsTests {

class FitTestSpacePoint {
 public:
  /// @brief Constructor for standard straw wires
  /// @param pos: Position of the wire
  /// @param driftR: Straw drift radius
  /// @param driftRUncert: Uncertainty on the drift radius uncertainty
  /// @param twinUncert: Uncertainty on the measurement along the straw
  FitTestSpacePoint(const Vector3& pos, const double driftR,
                    const double driftRUncert,
                    const std::optional<double> twinUncert = std::nullopt)
      : m_position{pos},
        m_driftR{driftR},
        m_measLoc0{twinUncert != std::nullopt} {
    using enum ResidualIdx;
    m_covariance[toUnderlying(bending)] = Acts::square(driftRUncert);
    m_covariance[toUnderlying(nonBending)] =
        Acts::square(twinUncert.value_or(0.));
  }
  /// @brief Constructor for rotated straw wires
  /// @param pos: Position of the wire
  /// @param wire: Orientation of the wire
  /// @param driftR: Drift radius of the measurement
  /// @param driftRUncert: Associated uncertainty on the measurement
  /// @param twinUncert: Uncertainty on the measurement along the straw
  FitTestSpacePoint(const Vector3& pos, const Vector3& wire,
                    const double driftR, const double driftRUncert,
                    const std::optional<double> twinUncert = std::nullopt)
      : m_position{pos},
        m_sensorDir{wire},
        m_driftR{driftR},
        m_measLoc0{twinUncert != std::nullopt} {
    using enum ResidualIdx;
    m_covariance[toUnderlying(bending)] = Acts::square(driftRUncert);
    m_covariance[toUnderlying(nonBending)] =
        Acts::square(twinUncert.value_or(0.));
  }

  /// @brief Constructor for spatial strip measurements
  /// @param stripPos: Position of the strip
  /// @param stripDir: Direction along the strip
  /// @param toNext: Vector pointing to the next strip inside the plane
  /// @param uncertLoc0: Uncertainty of the measurement in the non bending direction
  /// @param uncertLoc1: Uncertainty of the measurement in the bending direction
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
  /// @brief Constructor for strip measurements with time
  FitTestSpacePoint(const Vector3& stripPos, const Vector3& stripDir,
                    const Vector3& toNext, const double time,
                    const std::array<double, 3>& cov)
      : m_position{stripPos},
        m_sensorDir{stripDir},
        m_toNextSen{toNext},
        m_time{time},
        m_covariance{cov} {
    m_measLoc0 = m_covariance[toUnderlying(ResidualIdx::nonBending)] > 0.;
    m_measLoc1 = m_covariance[toUnderlying(ResidualIdx::bending)] > 0.;
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
  /// @brief Updates the time of the space point
  void updateTime(const double newTime) { m_time = newTime; }
  /// @brief Updates the status of the space point
  void updateStatus(const bool newStatus) { m_isGood = newStatus; }
  /// @brief Check if the measurement is valid after calibration
  bool isGood() const { return m_isGood; }
  /// @brief Returns the layer index of the space point
  std::size_t layer() const { return m_layer.value_or(0ul); }
  /// @brief Sets the layer number of the space point
  void setLayer(const std::size_t lay) {
    if (!m_layer) {
      m_layer = lay;
    }
  }

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
  bool m_isGood{true};
  std::optional<std::size_t> m_layer{};  // layer index starting from 0
};

static_assert(CompositeSpacePoint<FitTestSpacePoint>);

using Container_t = std::vector<std::shared_ptr<FitTestSpacePoint>>;
static_assert(CompositeSpacePointContainer<Container_t>);

class SpCalibrator {
 public:
  SpCalibrator() = default;

  explicit SpCalibrator(const Acts::Transform3& localToGlobal)
      : m_localToGlobal{localToGlobal} {}

  /// @brief Normalize input variable within a given range
  static double normalize(const double value, const double min,
                          const double max) {
    return 2. * (std::clamp(value, min, max) - min) / (max - min) - 1.;
  }
  /// @brief Drift time validity limits for the parametrized r(t) relation
  static constexpr double s_minDriftTime = 0._ns;
  static constexpr double s_maxDriftTime = 741.324 * 1._ns;
  /// @brief Normalize input drift time for r(t) relation
  static double normDriftTime(const double t) {
    return normalize(t, s_minDriftTime, s_maxDriftTime);
  }
  /// @brief Multiplicative factor arising from the derivative of the normalization
  static constexpr double s_timeNormFactor =
      2.0 / (s_maxDriftTime - s_minDriftTime);
  /// @brief Coefficients for the parametrized r(t) relation
  static constexpr std::array<double, 9> s_TtoRcoeffs = {
      9.0102,   6.7419,   -1.5829,   0.56475, -0.19245,
      0.019055, 0.030006, -0.034591, 0.023517};
  /// @brief Parametrized ATLAS r(t) relation
  static double driftRadius(const double t) {
    const double x = normDriftTime(t);
    double r{0.0};
    for (std::size_t n = 0; n < s_TtoRcoeffs.size(); ++n) {
      r += s_TtoRcoeffs[n] * Acts::detail::chebychevPolyTn(x, n);
    }
    return r;
  }
  /// @brief First derivative of r(t) relation
  static double driftVelocity(const double t) {
    const double x = normDriftTime(t);
    double v{0.0};
    for (std::size_t n = 1; n < s_TtoRcoeffs.size(); ++n) {
      v += s_TtoRcoeffs[n] * Acts::detail::chebychevPolyTn(x, n, 1u);
    }
    return v * s_timeNormFactor;
  }
  /// @brief Second derivative of r(t) relation
  static double driftAcceleration(const double t) {
    const double x = normDriftTime(t);
    double a{0.0};
    for (std::size_t n = 2; n < s_TtoRcoeffs.size(); ++n) {
      a += s_TtoRcoeffs[n] * Acts::detail::chebychevPolyTn(x, n, 2u);
    }
    return a * Acts::square(s_timeNormFactor);
  }
  /// @brief Drift radius validity limits for the parametrized t(r) relation
  static constexpr double s_minDriftRadius = 0.064502;
  static constexpr double s_maxDriftRadius = 14.566;
  /// @brief Normalize input drift radius for t(r) relation
  static double normDriftRadius(const double r) {
    return normalize(r, s_minDriftRadius, s_maxDriftRadius);
  }
  /// @brief Coefficients for the parametrized t(r) relation
  static constexpr std::array<double, 5> s_RtoTcoeffs = {
      256.476, 349.056, 118.349, 18.748, -6.4142};
  /// @brief Parametrized t(r) relation
  static double driftTime(const double r) {
    const double x = normDriftRadius(r);
    double t{0.0};
    for (std::size_t n = 0; n < s_RtoTcoeffs.size(); ++n) {
      t += s_RtoTcoeffs[n] * Acts::detail::legendrePoly(x, n);
    }
    return t * 1._ns;
  }
  /// @brief Coefficients for the uncertanty on the drift radius
  static constexpr std::array<double, 4> s_driftRUncertCoeffs{
      0.10826, -0.07182, 0.037597, -0.011712};
  /// @brief Compute the drift radius uncertanty
  static double driftUncert(const double r) {
    const double x = normDriftRadius(r);
    double s{0.0};
    for (std::size_t n = 0; n < s_driftRUncertCoeffs.size(); ++n) {
      s += s_driftRUncertCoeffs[n] * Acts::detail::chebychevPolyTn(x, n);
    }
    return s;
  }
  /// @brief Provide a fast estimate of the time of flight of the particle. Used in the Fast Fitter.
  /// @param measurement: measurement. It should be a straw measurement
  double fastToF(const FitTestSpacePoint& measurement) const {
    return (m_localToGlobal * measurement.localPosition()).norm() /
           PhysicalConstants::c;
  }
  /// @brief Provide the calibrated drift radius given the straw measurement and time offset
  ///        Needed for the Fast Fitter.
  /// @param ctx: Calibration context (Needed by concept interface)
  /// @param measurement: measurement. It should be a straw measurement
  /// @param timeOffSet: Offset in the time of arrival
  double driftRadius(const Acts::CalibrationContext& /*ctx*/,
                     const FitTestSpacePoint& measurement,
                     const double timeOffSet) const {
    if (!measurement.isStraw() || !measurement.isGood()) {
      return 0.;
    }
    return driftRadius(measurement.time() - timeOffSet - fastToF(measurement));
  }
  /// @brief Provide the drift velocity given the straw measurent and time offset
  ///        Needed for the Fast Fitter.
  /// @param ctx: Calibration context (Needed by concept interface)
  /// @param measurement: measurement
  /// @param timeOffSet: Offset in the time of arrival
  double driftVelocity(const Acts::CalibrationContext& /*ctx*/,
                       const FitTestSpacePoint& measurement,
                       const double timeOffSet) const {
    if (!measurement.isStraw() || !measurement.isGood()) {
      return 0.;
    }
    return driftVelocity(measurement.time() - timeOffSet -
                         fastToF(measurement));
  }
  /// @brief Provide the drift acceleration given the straw measurent and time offset
  ///        Needed for the Fast Fitter.
  /// @param ctx: Calibration context (Needed by concept interface)
  /// @param measurement: measurement
  /// @param timeOffSet: Offset in the time of arrival
  double driftAcceleration(const Acts::CalibrationContext& /*ctx*/,
                           const FitTestSpacePoint& measurement,
                           const double timeOffSet) const {
    if (!measurement.isStraw() || !measurement.isGood()) {
      return 0.;
    }
    return driftAcceleration(measurement.time() - timeOffSet -
                             fastToF(measurement));
  }
  /// @brief Compute the distance of the point of closest approach of a straw measurement
  /// @param ctx: Calibration context (Needed by concept interface)
  /// @param trackPos: Position of the track at z=0.
  /// @param trackDir: Direction of the track in the local frame
  /// @param measurement: Measurement
  double closestApproachDist(const Vector3& trackPos, const Vector3& trackDir,
                             const FitTestSpacePoint& measurement) const {
    const Vector3 closePointOnLine =
        lineIntersect(measurement.localPosition(),
                      measurement.sensorDirection(), trackPos, trackDir)
            .position();
    return (m_localToGlobal * closePointOnLine).norm();
  }
  std::shared_ptr<FitTestSpacePoint> calibrate(
      const Acts::CalibrationContext& /*cctx*/, const Vector3& trackPos,
      const Vector3& trackDir, const double timeOffSet,
      const FitTestSpacePoint& sp) const {
    auto copySp = std::make_shared<FitTestSpacePoint>(sp);
    if (!copySp->measuresLoc0() || !copySp->measuresLoc1()) {
      /// Estimate the best position along the sensor
      auto bestPos = lineIntersect(trackPos, trackDir, copySp->localPosition(),
                                   copySp->sensorDirection());
      copySp->updatePosition(bestPos.position());
    }
    if (copySp->isStraw()) {
      const double driftTime{copySp->time() - timeOffSet -
                             closestApproachDist(trackPos, trackDir, sp) /
                                 PhysicalConstants::c};
      if (driftTime > s_minDriftTime && driftTime < s_maxDriftTime) {
        copySp->updateStatus(true);
        copySp->updateDriftR(driftRadius(driftTime));
      } else {
        copySp->updateStatus(false);
        copySp->updateDriftR(0.);
      }
    }
    return copySp;
  }
  /// @brief Calibrate a set of straw measurements using the best known estimate on a straight line track
  /// @param ctx: Calibration context (Needed by concept interface)
  /// @param trackPos: Position of the track at z=0.
  /// @param trackDir: Direction of the track in the local frame
  /// @param timeOffSet: Offset in the time of arrival (To be implemented)
  /// @param uncalibCont: Uncalibrated composite space point container
  Container_t calibrate(const Acts::CalibrationContext& ctx,
                        const Vector3& trackPos, const Vector3& trackDir,
                        const double timeOffSet,
                        const Container_t& uncalibCont) const {
    Container_t calibMeas{};
    std::ranges::transform(
        uncalibCont, std::back_inserter(calibMeas), [&](const auto& calibMe) {
          return calibrate(ctx, trackPos, trackDir, timeOffSet, *calibMe);
        });
    return calibMeas;
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
  /// @brief Provide the drift velocity given the straw measurent after being calibrated
  /// @param ctx: Calibration context (Needed by concept interface)
  /// @param measurement: measurement
  static double driftVelocity(const Acts::CalibrationContext& /*ctx*/,
                              const FitTestSpacePoint& measurement) {
    if (!measurement.isStraw() || !measurement.isGood()) {
      return 0.;
    }
    return driftVelocity(driftTime(Acts::abs(measurement.driftRadius())));
  }
  /// @brief Provide the drift acceleration given the straw measurent after being calibrated
  /// @param ctx: Calibration context (Needed by concept interface)
  /// @param measurement: measurement
  static double driftAcceleration(const Acts::CalibrationContext& /*ctx*/,
                                  const FitTestSpacePoint& measurement) {
    if (!measurement.isStraw() || !measurement.isGood()) {
      return 0.;
    }
    return driftAcceleration(driftTime(Acts::abs(measurement.driftRadius())));
  }

 private:
  const Acts::Transform3 m_localToGlobal{Acts::Transform3::Identity()};
};
/// Ensure that the Test space point calibrator satisfies the calibrator concept
static_assert(
    CompositeSpacePointCalibrator<SpCalibrator, Container_t, Container_t>);
static_assert(
    CompositeSpacePointFastCalibrator<SpCalibrator, FitTestSpacePoint>);

/// @brief Split the composite space point container into straw and strip measurements
class SpSorter {
 protected:
  /// @brief Protected constructor to instantiate the SpSorter.
  /// @param hits: List of space points to sort per layer
  /// @param calibrator: space point calibrator to recalibrate the hits
  SpSorter(const Container_t& hits, const SpCalibrator* calibrator)
      : m_calibrator{calibrator} {
    for (const auto& spPtr : hits) {
      auto& pushMe{spPtr->isStraw() ? m_straws : m_strips};
      if (spPtr->layer() >= pushMe.size()) {
        pushMe.resize(spPtr->layer() + 1);
      }
      pushMe[spPtr->layer()].push_back(spPtr);
    }
  }

 public:
  /// @brief Return the sorted straw hits
  const std::vector<Container_t>& strawHits() const { return m_straws; }
  /// @brief Returns the sorted strip hits
  const std::vector<Container_t>& stripHits() const { return m_strips; }
  /// @brief Returns whether the candidate space point is a good hit or not
  /// @param testPoint: Hit to test
  bool goodCandidate(const FitTestSpacePoint& testPoint) const {
    return testPoint.isGood();
  }
  /// @brief Calculates the pull of the space point w.r.t. to the
  ///        candidate seed line. To improve the pull's precision
  ///        the function may call the calibrator in the backend
  /// @param cctx: Reference to the calibration context to pipe
  ///              the hook for conditions access to the caller
  /// @param pos: Position of the cancidate seed line
  /// @param dir: Direction of the candidate seed line
  /// @param t0: Offse in the time of arrival of the particle
  /// @param testSp: Reference to the straw space point to test
  double candidateChi2(const CalibrationContext& /* cctx*/, const Vector3& pos,
                       const Vector3& dir, const double /*t0*/,
                       const FitTestSpacePoint& testSp) const {
    return CompSpacePointAuxiliaries::chi2Term(pos, dir, testSp);
  }
  /// @brief Creates a new empty container
  Container_t newContainer(const CalibrationContext& /*cctx*/) const {
    return Container_t{};
  }
  /// @brief Appends the space point to the container and calibrates it according to
  ///        the seed parameters
  ///
  void append(const CalibrationContext& cctx, const Vector3& pos,
              const Vector3& dir, const double t0,
              const FitTestSpacePoint& testSp,
              Container_t& outContainer) const {
    outContainer.push_back(m_calibrator->calibrate(cctx, pos, dir, t0, testSp));
  }

  bool stopSeeding(const std::size_t lowerLayer,
                   const std::size_t upperLayer) const {
    return lowerLayer >= upperLayer;
  }
  double strawRadius(const FitTestSpacePoint& /*testSp*/) const {
    return 15._mm;
  }

 private:
  const SpCalibrator* m_calibrator{nullptr};
  std::vector<Container_t> m_straws{};
  std::vector<Container_t> m_strips{};
};

static_assert(CompositeSpacePointSorter<SpSorter, Container_t>);

/// @brief Generates a random straight line
/// @param engine Random number sequence to draw the parameters from
/// @return A Line object instantiated with the generated parameters
Line_t generateLine(RandomEngine& engine, const Logger& logger) {
  using ParIndex = Line_t::ParIndex;
  Line_t::ParamVector linePars{};
  linePars[toUnderlying(ParIndex::phi)] =
      uniform{-120_degree, 120_degree}(engine);
  linePars[toUnderlying(ParIndex::x0)] = uniform{-500., 500.}(engine);
  linePars[toUnderlying(ParIndex::y0)] = uniform{-500., 500.}(engine);
  linePars[toUnderlying(ParIndex::theta)] =
      uniform{5_degree, 175_degree}(engine);
  if ((Acts::abs(linePars[toUnderlying(ParIndex::theta)] - 90._degree) <
       10._degree) ||
      (Acts::abs(linePars[toUnderlying(ParIndex::phi)]) < 15._degree)) {
    return generateLine(engine, logger);
  }
  Line_t line{linePars};

  ACTS_DEBUG(
      "\n\n\n"
      << __func__ << "() " << __LINE__ << " - Generated parameters "
      << std::format("theta: {:.2f}, phi: {:.2f}, y0: {:.1f}, x0: {:.1f}",
                     linePars[toUnderlying(ParIndex::theta)] / 1._degree,
                     linePars[toUnderlying(ParIndex::phi)] / 1._degree,
                     linePars[toUnderlying(ParIndex::y0)],
                     linePars[toUnderlying(ParIndex::x0)])
      << " --> " << toString(line.position()) << " + "
      << toString(line.direction()));

  return line;
}

class MeasurementGenerator {
 public:
  struct Config {
    /// @brief Create straw measurements
    bool createStraws{true};
    /// @brief Smear the straw radius
    bool smearRadius{true};
    /// @brief Straw measurement measures the coordinate along the wire
    bool twinStraw{false};
    /// @brief Resolution of the coordinate along the wire measurement
    double twinStrawReso{5._cm};
    /// @brief Create strip measurements
    bool createStrips{true};
    /// @brief Smear the strips around the pitch
    bool smearStrips{true};
    /// @brief Alternatively, discretize the strips onto a strip plane
    bool discretizeStrips{false};
    /// @brief Combine the two strip measurements to a single space point
    bool combineSpacePoints{false};
    /// @brief Pitch between two loc0 strips
    double stripPitchLoc0{4._cm};
    /// @brief Pitch between two loc1 strips
    double stripPitchLoc1{3._cm};
    /// @brief Direction of the strip if it measures loc1
    std::vector<Vector3> stripDirLoc1 =
        std::vector<Vector3>(8, Vector3::UnitX());
    /// @brief Direction of the strip if it measures loc0
    std::vector<Vector3> stripDirLoc0 =
        std::vector<Vector3>(8, Vector3::UnitY());
    /// @brief Experiment specific calibration context
    Acts::CalibrationContext calibContext{};
    /// @brief Local to global transform
    Acts::Transform3 localToGlobal{Acts::Transform3::Identity()};
  };
  /// @brief Extrapolate the straight line track through the straw layers to
  ///        evaluate which tubes were crossed by the track. The straw layers
  ///        are staggered in the z-direction and each layer expands in y. To
  ///        estimate, which tubes were crossed, the track is extrapolated to
  ///        the z-plane below & above the straw wires. Optionally, the true
  ///        drift-radius can be smeared assuming a Gaussian with a drift-radius
  ///        dependent uncertainty.
  /// @param line: The track to extrapolate
  /// @param engine: Random number generator to smear the drift radius
  /// @param smearRadius: If true, the drift radius is smeared with a Gaussian
  /// @param createStrips: If true the strip measurements are created
  static Container_t spawn(const Line_t& line, const double t0,
                           RandomEngine& engine, const Config& genCfg,
                           const Logger& logger) {
    /// Direction vector to go a positive step in the tube honeycomb grid
    const Vector3 posStaggering{0., std::cos(60._degree), std::sin(60._degree)};
    /// Direction vector to go a negative step in the tube honeycomb grid
    const Vector3 negStaggering{0., -std::cos(60._degree),
                                std::sin(60._degree)};
    /// Number of tube layers per multilayer
    constexpr std::size_t nLayersPerMl = 8;
    /// Number of overall tubelayers
    constexpr std::size_t nTubeLayers = nLayersPerMl * 2;
    /// Position in z of the first tube layer
    constexpr double chamberDistance = 5._m;
    /// Radius of each straw
    constexpr double tubeRadius = 15._mm;
    /// Distance between the first <nLayersPerMl> layers and the second pack
    constexpr double tubeLayerDist = 1.2_m;
    /// Distance between the tube multilayers and to the first strip layer
    constexpr double tubeStripDist = 30._cm;
    /// Distance between two strip layers
    constexpr double stripLayDist = 0.5_cm;

    std::array<Vector3, nTubeLayers> tubePositions{
        filledArray<Vector3, nTubeLayers>(chamberDistance * Vector3::UnitZ())};
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
    if (genCfg.createStraws) {
      const SpCalibrator calibrator{genCfg.localToGlobal};
      for (const auto& [i_layer, stag] : enumerate(tubePositions)) {
        auto planeExtpLow =
            intersectPlane(line.position(), line.direction(), Vector3::UnitZ(),
                           stag.z() - tubeRadius);
        auto planeExtpHigh =
            intersectPlane(line.position(), line.direction(), Vector3::UnitZ(),
                           stag.z() + tubeRadius);
        ACTS_DEBUG("spawn() - Extrapolated to plane "
                   << toString(planeExtpLow.position()) << " "
                   << toString(planeExtpHigh.position()));

        const auto dToFirstLow = static_cast<int>(std::ceil(
            (planeExtpLow.position().y() - stag.y()) / (2. * tubeRadius)));
        const auto dToFirstHigh = static_cast<int>(std::ceil(
            (planeExtpHigh.position().y() - stag.y()) / (2. * tubeRadius)));
        /// Does the track go from left to right or right to left?
        const int dT = copySign(1, dToFirstHigh - dToFirstLow);
        /// Loop over the candidate tubes and check each one whether the track
        /// actually crossed them. Then generate the circle and optionally smear
        /// the radius
        for (int tN = dToFirstLow - dT; tN != dToFirstHigh + 2 * dT; tN += dT) {
          Vector3 tube = stag + 2. * tN * tubeRadius * Vector3::UnitY();
          const double rad = Acts::abs(signedDistance(
              tube, Vector3::UnitX(), line.position(), line.direction()));
          if (rad > tubeRadius) {
            continue;
          }

          const double smearedR =
              genCfg.smearRadius
                  ? Acts::abs(
                        normal_t{rad, SpCalibrator::driftUncert(rad)}(engine))
                  : rad;
          if (smearedR > tubeRadius) {
            continue;
          }
          ///
          if (genCfg.twinStraw) {
            const auto iSectWire = lineIntersect<3>(
                line.position(), line.direction(), tube, Vector3::UnitX());
            tube = iSectWire.position();
            tube[eX] = normal_t{tube[eX], genCfg.twinStrawReso}(engine);
          }
          ACTS_DEBUG("spawn() - Tube position: " << toString(tube)
                                                 << ", radius: " << rad);

          auto twinUncert =
              genCfg.twinStraw
                  ? std::make_optional<double>(genCfg.twinStrawReso)
                  : std::nullopt;
          auto& sp =
              measurements.emplace_back(std::make_unique<FitTestSpacePoint>(
                  tube, smearedR, SpCalibrator::driftUncert(smearedR),
                  twinUncert));
          sp->setLayer(i_layer);
          const double driftTime{SpCalibrator::driftTime(sp->driftRadius())};
          const double tof{calibrator.closestApproachDist(
                               line.position(), line.direction(), *sp) /
                           PhysicalConstants::c};
          sp->updateTime(driftTime + t0 + tof);
          ACTS_VERBOSE("spawn() - Created "
                       << toString(*sp) << ", t: " << sp->time() / 1._ns
                       << " ns"
                       << ", tDrift: " << driftTime / 1._ns << " ns"
                       << ", ToF: " << tof / 1._ns << " ns.");
        }
      }
    }
    if (genCfg.createStrips) {
      ///@brief Helper function to discretize the extrapolated position on the strip plane to actual strips.
      ///       The function returns the local coordinate in the picked plane
      ///       projection
      /// @param extPos: Extrapolated position of the straight line in the plane
      /// @param loc1: Boolean to fetch either the loc1 projection or the loc0 projection
      auto discretize = [&genCfg, &engine](const Vector3& extPos,
                                           const std::size_t lay,
                                           const bool loc1) -> Vector3 {
        const double pitch =
            loc1 ? genCfg.stripPitchLoc1 : genCfg.stripPitchLoc0;
        assert(pitch > 0.);

        const Vector3& stripDir =
            loc1 ? genCfg.stripDirLoc1.at(lay) : genCfg.stripDirLoc0.at(lay);
        const Vector3 stripNormal = stripDir.cross(Vector3::UnitZ());
        const double dist = stripNormal.dot(extPos);
        if (genCfg.smearStrips) {
          return normal_t{dist, pitch / std::sqrt(12)}(engine)*stripNormal;
        }
        if (genCfg.discretizeStrips) {
          return pitch *
                 static_cast<int>(std::ceil((dist - 0.5 * pitch) / pitch)) *
                 stripNormal;
        }
        return dist * stripNormal;
      };
      /// Calculate the strip measurement's covariance
      const double stripCovLoc0 =
          Acts::square(genCfg.stripPitchLoc0) / std::sqrt(12.);
      const double stripCovLoc1 =
          Acts::square(genCfg.stripPitchLoc1) / std::sqrt(12.);

      for (std::size_t sL = 0; sL < std::max(genCfg.stripDirLoc0.size(),
                                             genCfg.stripDirLoc1.size());
           ++sL) {
        const double distInZ = tubeStripDist + sL * stripLayDist;
        const double planeLow = tubePositions.front().z() - distInZ;
        const double planeHigh = tubePositions.back().z() + distInZ;

        for (const double plane : {planeLow, planeHigh}) {
          const auto extp = intersectPlane(line.position(), line.direction(),
                                           Vector3::UnitZ(), plane);
          ACTS_VERBOSE("spawn() - Propagated line to "
                       << toString(extp.position()) << ".");

          if (genCfg.combineSpacePoints &&
              sL < std::min(genCfg.stripDirLoc0.size(),
                            genCfg.stripDirLoc1.size())) {
            const Vector3 extpPos{discretize(extp.position(), sL, false) +
                                  discretize(extp.position(), sL, true) +
                                  plane * Vector3::UnitZ()};
            measurements.emplace_back(std::make_unique<FitTestSpacePoint>(
                extpPos, genCfg.stripDirLoc0.at(sL), genCfg.stripDirLoc1.at(sL),
                stripCovLoc0, stripCovLoc1));
            measurements.back()->setLayer(sL);
          } else {
            if (sL < genCfg.stripDirLoc0.size()) {
              const Vector3 extpPos{discretize(extp.position(), sL, false) +
                                    plane * Vector3::UnitZ()};
              measurements.emplace_back(std::make_unique<FitTestSpacePoint>(
                  extpPos, genCfg.stripDirLoc0.at(sL),
                  genCfg.stripDirLoc0.at(sL).cross(Vector3::UnitZ()),
                  stripCovLoc0, 0.));
              auto& nM{*measurements.back()};
              nM.setLayer(sL);
              ACTS_VERBOSE("spawn() - Created loc0 strip @" << toString(nM)
                                                            << ".");
            }
            if (sL < genCfg.stripDirLoc1.size()) {
              const Vector3 extpPos{discretize(extp.position(), sL, true) +
                                    plane * Vector3::UnitZ()};
              measurements.emplace_back(std::make_unique<FitTestSpacePoint>(
                  extpPos, genCfg.stripDirLoc1.at(sL),
                  genCfg.stripDirLoc1.at(sL).cross(Vector3::UnitZ()), 0.,
                  stripCovLoc1));
              auto& nM{*measurements.back()};
              nM.setLayer(sL);
              ACTS_VERBOSE("spawn() - Created loc1 strip @" << toString(nM)
                                                            << ".");
            }
          }
        }
      }
    }
    ACTS_DEBUG("Track hit in total " << measurements.size() << " sensors.");
    std::ranges::sort(measurements, [&line](const auto& spA, const auto& spB) {
      return line.direction().dot(spA->localPosition() - line.position()) <
             line.direction().dot(spB->localPosition() - line.position());
    });
    return measurements;
  }
};

/// @brief Construct the start parameters from the hit container && the true trajectory
/// @param line: True trajectory to pick-up the correct left/right ambiguity for the straw seed hits
/// @param hits: List of measurements to be used for fitting
ParamVec_t startParameters(const Line_t& line, const Container_t& hits) {
  ParamVec_t pars{};

  double tanAlpha{0.};
  double tanBeta{0.};
  /// Setup the seed parameters in x0 && phi
  auto firstPhi = std::ranges::find_if(
      hits, [](const auto& sp) { return sp->measuresLoc0(); });
  auto lastPhi =
      std::ranges::find_if(std::ranges::reverse_view(hits),
                           [](const auto& sp) { return sp->measuresLoc0(); });

  if (firstPhi != hits.end() && lastPhi != hits.rend()) {
    const Vector3 firstToLastPhi =
        (**lastPhi).localPosition() - (**firstPhi).localPosition();
    tanAlpha = firstToLastPhi.x() / firstToLastPhi.z();
    /// -> x = tanPhi * z + x_{0} ->
    pars[toUnderlying(FitParIndex::x0)] =
        (**lastPhi).localPosition().x() -
        (**lastPhi).localPosition().z() * tanAlpha;
  }
  /// Setup the seed parameters in y0 && theta
  auto firstTube =
      std::ranges::find_if(hits, [](const auto& sp) { return sp->isStraw(); });
  auto lastTube =
      std::ranges::find_if(std::ranges::reverse_view(hits),
                           [](const auto& sp) { return sp->isStraw(); });

  if (firstTube != hits.end() && lastTube != hits.rend()) {
    const int signFirst =
        CompSpacePointAuxiliaries::strawSign(line, **firstTube);
    const int signLast = CompSpacePointAuxiliaries::strawSign(line, **lastTube);

    auto seedPars = CompositeSpacePointLineSeeder::constructTangentLine(
        **lastTube, **firstTube,
        CompositeSpacePointLineSeeder::encodeAmbiguity(signLast, signFirst));
    tanBeta = std::tan(seedPars.theta);
    pars[toUnderlying(FitParIndex::y0)] = seedPars.y0;
  } else {
    auto firstEta = std::ranges::find_if(hits, [](const auto& sp) {
      return !sp->isStraw() && sp->measuresLoc1();
    });
    auto lastEta = std::ranges::find_if(
        std::ranges::reverse_view(hits),
        [](const auto& sp) { return !sp->isStraw() && sp->measuresLoc1(); });

    if (firstEta != hits.end() && lastEta != hits.rend()) {
      const Vector3 firstToLastEta =
          (**lastEta).localPosition() - (**firstEta).localPosition();
      tanBeta = firstToLastEta.y() / firstToLastEta.z();
      /// -> y = tanTheta * z + y_{0} ->
      pars[toUnderlying(FitParIndex::y0)] =
          (**lastEta).localPosition().y() -
          (**lastEta).localPosition().z() * tanBeta;
    }
  }
  const Vector3 seedDir = makeDirectionFromAxisTangents(tanAlpha, tanBeta);

  pars[toUnderlying(FitParIndex::theta)] = theta(seedDir);
  pars[toUnderlying(FitParIndex::phi)] = phi(seedDir);
  return pars;
}

bool isGoodHit(const FitTestSpacePoint& sp) {
  return !sp.isStraw() || sp.isGood();
};

}  // namespace ActsTests

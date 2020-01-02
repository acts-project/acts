// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/EventData/TrackStateSorters.hpp"
#include "Acts/Fitter/KalmanFitterError.hpp"
#include "Acts/Fitter/detail/VoidKalmanComponents.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <functional>
#include <map>
#include <memory>

namespace Acts {

/// @brief Options struct how the Fitter is called
///
/// It contains the context of the fitter call and the optional
/// surface where to express the fit result
///
/// @note the context objects must be provided
struct TrackFinderOptions {
  /// Deleted default constructor
  TrackFinderOptions() = delete;

  /// PropagatorOptions with context
  ///
  /// @param gctx The goemetry context for this fit
  /// @param mctx The magnetic context for this fit
  /// @param cctx The calibration context for this fit
  /// @param rSurface The reference surface for the fit to be expressed at
  TrackFinderOptions(std::reference_wrapper<const GeometryContext> gctx,
                     std::reference_wrapper<const MagneticFieldContext> mctx,
                     std::reference_wrapper<const CalibrationContext> cctx,
                     const Surface* rSurface = nullptr, bool mScattering = true,
                     bool eLoss = true)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        referenceSurface(rSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss) {}

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  /// The reference Surface
  const Surface* referenceSurface = nullptr;

  /// Whether to consider multiple scattering.
  bool multipleScattering = true;

  /// Whether to consider energy loss.
  bool energyLoss = true;
};

template <typename source_link_t>
struct TrackFinderResult {
  // Fitted states that the actor has handled.
  MultiTrajectory<source_link_t> fittedStates;

  // This is the index of the 'tip' of the track stored in multitrajectory.
  // Since this KF only stores one trajectory, it is unambiguous.
  // SIZE_MAX is the start of a trajectory.
  size_t trackTip = SIZE_MAX;

  // Counter for states with measurements
  size_t measurementStates = 0;

  // Counter for handled states
  size_t processedStates = 0;

  // Indicator if initialization has been performed.
  bool initialized = false;

  Result<void> result{Result<void>::success()};
};

/// @brief Track finder implementation of Acts as a plugin
///
/// to the Propgator
///
/// @tparam propagator_t Type of the propagation class
/// @tparam updater_t Type of the kalman updater class
/// @tparam smoother_t Type of the kalman smoother class
/// @tparam calibrator_t Type of the calibrator class
/// @tparam input_converter_t Type of the input converter class
/// @tparam output_converter_t Type of the output converter class
///
/// The track finder contains an Actor and a Sequencer sub-class.
/// The Sequencer has to be part of the Navigator of the Propagator
/// in order to initialize and provide the measurement surfaces.
///
/// The Actor is part of the Propagation call and does the Kalman update
/// and eventually the smoothing.  Updater, Smoother and Calibrator are
/// given to the Actor for further use:
/// - The Updater is the implemented kalman updater formalism, it
///   runs via a visitor pattern through the measurements.
/// - The Smoother is called at the end of the forward fit by the Actor.
/// - The Calibrator is a dedicated calibration algorithm that allows
///   to calibrate measurements using track information, this could be
///    e.g. sagging for wires, module deformations, etc.
///
/// Measurements are not required to be ordered for the track finder,
/// measurement ordering needs to be figured out by the navigation of
/// the propagator.
///
/// The Input converter is a converter that transforms the input
/// measurement/track/segments into a set of FittableMeasurements
///
/// The Output converter is a converter that transforms the
/// set of track states into a given track/track particle class
///
/// The void components are provided mainly for unit testing.
template <typename propagator_t, typename updater_t = VoidKalmanUpdater,
          typename smoother_t = VoidKalmanSmoother,
          typename calibrator_t = VoidMeasurementCalibrator,
          typename input_converter_t = VoidKalmanComponents,
          typename output_converter_t = VoidKalmanComponents>
class TrackFinder {
 public:
  /// Shorthand definition
  using MeasurementSurfaces = std::multimap<const Layer*, const Surface*>;

  /// Default constructor is deleted
  TrackFinder() = delete;

  /// Constructor from arguments
  TrackFinder(propagator_t pPropagator,
              std::unique_ptr<const Logger> logger =
                  getDefaultLogger("TrackFinder", Logging::INFO),
              input_converter_t pInputCnv = input_converter_t(),
              output_converter_t pOutputCnv = output_converter_t())
      : m_propagator(std::move(pPropagator)),
        m_inputConverter(std::move(pInputCnv)),
        m_outputConverter(std::move(pOutputCnv)),
        m_logger(logger.release()) {}

 private:
  /// The propgator for the transport and material update
  propagator_t m_propagator;

  /// The input converter to Fittable measurements
  input_converter_t m_inputConverter;

  /// The output converter into a given format
  output_converter_t m_outputConverter;

  /// Logger getter to support macros
  const Logger& logger() const { return *m_logger; }

  /// Owned logging instance
  std::shared_ptr<const Logger> m_logger;

  /// The navigator type
  using KalmanNavigator = typename decltype(m_propagator)::Navigator;

  /// The navigator has DirectNavigator type or not
  static constexpr bool isDirectNavigator =
      std::is_same<KalmanNavigator, DirectNavigator>::value;

  /// @brief Propagator Actor plugin for the TrackFinder
  ///
  /// @tparam source_link_t is an type fulfilling the @c SourceLinkConcept
  /// @tparam parameters_t The type of parameters used for "local" paremeters.
  ///
  /// The TrackFinderActor does not rely on the measurements to be
  /// sorted along the track.
  template <typename source_link_t, typename parameters_t>
  class Actor {
   public:
    using TrackStateType = TrackState<source_link_t, parameters_t>;

    /// Explicit constructor with updater and calibrator
    Actor(updater_t pUpdater = updater_t(), smoother_t pSmoother = smoother_t(),
          calibrator_t pCalibrator = calibrator_t())
        : m_updater(std::move(pUpdater)),
          m_smoother(std::move(pSmoother)),
          m_calibrator(std::move(pCalibrator)) {}

    /// Broadcast the result_type
    using result_type = TrackFinderResult<source_link_t>;

    /// The target surface
    const Surface* targetSurface = nullptr;

    /// Allows retrieving measurements for a surface
    std::multimap<const Surface*, source_link_t> inputMeasurements;

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// Whether to consider energy loss.
    bool energyLoss = true;

    /// @brief Track finder actor operation
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    void operator()(propagator_state_t& state, const stepper_t& stepper,
                    result_type& result) const {
      ACTS_VERBOSE("TrackFinder step");
      // Initialization:
      // - Only when track states are not set
      if (!result.initialized) {
        // -> Move the TrackState vector
        // -> Feed the KalmanSequencer with the measurements to be fitted
        ACTS_VERBOSE("Initializing");
        initialize(state, stepper, result);
        result.initialized = true;
      }

      // Update:
      // - Waiting for a current surface that has material
      // -> a trackState will be created on surface with material
      auto surface = state.navigation.currentSurface;
      if (surface and surface->surfaceMaterial() and
          not state.navigation.targetReached) {
        // Check if the surface is in the measurement map
        // -> Get the measurement / calibrate
        // -> Create the predicted state
        // -> Perform the kalman update
        // -> Check outlier behavior (@todo)
        // -> Fill strack state information & update stepper information
        ACTS_VERBOSE("Perform filter step");
        auto res = filter(surface, state, stepper, result);
        if (!res.ok()) {
          ACTS_ERROR("Error in filter: " << res.error());
          result.result = res.error();
        }
      }

      // Finalization:
      // When all track states have been handled or the navigation is breaked
      if (result.measurementStates == multimapKeyCount(inputMeasurements) or
          (result.measurementStates > 0 and state.navigation.navigationBreak)) {
        // Manually break the navigation here
        state.navigation.targetReached = true;
      }
    }

    /// @brief Track finder actor operation : initialize
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    void initialize(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    result_type& /*result*/) const {}

    /// @brief Track finder actor operation : update
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param result The mutable result state object
    template <typename propagator_state_t, typename stepper_t>
    Result<void> filter(const Surface* surface, propagator_state_t& state,
                        const stepper_t& stepper, result_type& result) const {
      // Try to find the surface in the measurement surfaces
      auto sourcelink_it = inputMeasurements.find(surface);
      if (sourcelink_it != inputMeasurements.end()) {
        // Screen output message
        ACTS_VERBOSE("Measurement surface " << surface->geoID()
                                            << " detected.");

        // Update state and stepper with pre material effects
        materialInteractor(surface, state, stepper, preUpdate);

        // Transport & bind the state to the current surface
        auto [boundParams, jacobian, pathLength] =
            stepper.boundState(state.stepping, *surface, true);

        // add a full TrackState entry multi trajectory
        // (this allocates storage for all components, we will set them later)
        result.trackTip = result.fittedStates.addTrackState(
            TrackStatePropMask::All, result.trackTip);

        // now get track state proxy back
        auto trackStateProxy =
            result.fittedStates.getTrackState(result.trackTip);

        // Fill the track state
        trackStateProxy.predicted() = boundParams.parameters();
        trackStateProxy.predictedCovariance() = *boundParams.covariance();
        trackStateProxy.jacobian() = jacobian;
        trackStateProxy.pathLength() = pathLength;

        // Get and set the type flags
        auto& typeFlags = trackStateProxy.typeFlags();
        typeFlags.set(TrackStateFlag::MaterialFlag);
        typeFlags.set(TrackStateFlag::MeasurementFlag);
        typeFlags.set(TrackStateFlag::ParameterFlag);

        // Create a tempoary track state proxy
        size_t tempTrackTip =
            result.fittedStates.addTrackState(TrackStatePropMask::All);

        // Get the temporary track state proxy
        auto tempTrackStateProxy =
            result.fittedStates.getTrackState(tempTrackTip);

        // Fill the temporary track state proxy
        tempTrackStateProxy.predicted() = boundParams.parameters();
        tempTrackStateProxy.predictedCovariance() = *boundParams.covariance();
        tempTrackStateProxy.jacobian() = jacobian;
        tempTrackStateProxy.pathLength() = pathLength;

        // Initialize the update result
        Result<void> updateRes = KalmanFitterError::UpdateFailed;

        // Get all the source links on this surface
        auto sourcelinks = inputMeasurements.equal_range(surface);

        // Loop over the source link to select candidate
        // with the minimum chi2
        double minChi2 = 999;
        for (auto it = sourcelinks.first; it != sourcelinks.second; ++it) {
          // assign the source link to the track state
          tempTrackStateProxy.uncalibrated() = it->second;

          // We have predicted parameters, so calibrate the uncalibrated input
          // measuerement
          std::visit(
              [&](const auto& calibrated) {
                tempTrackStateProxy.setCalibrated(calibrated);
              },
              m_calibrator(tempTrackStateProxy.uncalibrated(),
                           tempTrackStateProxy.predicted()));

          // If the update is successful, set covariance and
          auto tempUpdateRes = m_updater(state.geoContext, tempTrackStateProxy);
          if (tempUpdateRes.ok()) {
            if (tempTrackStateProxy.chi2() < minChi2) {
              updateRes = tempUpdateRes;
              minChi2 = tempTrackStateProxy.chi2();
              // Update the track state proxy data with the tempoary track state
              // proxy
              trackStateProxy.uncalibrated() = it->second;
              std::visit(
                  [&](const auto& calibrated) {
                    trackStateProxy.setCalibrated(calibrated);
                  },
                  m_calibrator(trackStateProxy.uncalibrated(),
                               trackStateProxy.predicted()));
              trackStateProxy.filtered() = tempTrackStateProxy.filtered();
              trackStateProxy.filteredCovariance() =
                  tempTrackStateProxy.filteredCovariance();
              trackStateProxy.chi2() = tempTrackStateProxy.chi2();
            }
          }
        }  // end of loop for all source links on this surface

        // If the update is successful, set covariance and
        if (!updateRes.ok()) {
          ACTS_ERROR("Update step failed: " << updateRes.error());
          return updateRes.error();
        } else {
          // Update the stepping state with filtered parameters
          ACTS_VERBOSE("Filtering step successful, updated parameters are : \n"
                       << trackStateProxy.filtered().transpose());
          // update stepping state using filtered parameters after kalman update
          // We need to (re-)construct a BoundParameters instance here, which is
          // a bit awkward.
          stepper.update(state.stepping, trackStateProxy.filteredParameters(
                                             state.options.geoContext));

          // Update state and stepper with post material effects
          materialInteractor(surface, state, stepper, postUpdate);
        }
        // We count the state with measurement
        ++result.measurementStates;
        // We count the processed state
        ++result.processedStates;
      } else {
        // add a non-measurement TrackState entry multi trajectory
        // (this allocates storage for components except measurements, we will
        // set them later)
        result.trackTip = result.fittedStates.addTrackState(
            ~(TrackStatePropMask::Uncalibrated |
              TrackStatePropMask::Calibrated),
            result.trackTip);

        // now get track state proxy back
        auto trackStateProxy =
            result.fittedStates.getTrackState(result.trackTip);

        // Set the surface
        trackStateProxy.setReferenceSurface(surface->getSharedPtr());

        // Set the track state flags
        auto& typeFlags = trackStateProxy.typeFlags();
        typeFlags.set(TrackStateFlag::MaterialFlag);
        typeFlags.set(TrackStateFlag::ParameterFlag);

        if (surface->associatedDetectorElement() != nullptr) {
          ACTS_VERBOSE("Detected hole on " << surface->geoID());
          // If the surface is sensitive, set the hole type flag
          typeFlags.set(TrackStateFlag::HoleFlag);

          // Transport & bind the state to the current surface
          auto [boundParams, jacobian, pathLength] =
              stepper.boundState(state.stepping, *surface, true);

          // Fill the track state
          trackStateProxy.predicted() = boundParams.parameters();
          trackStateProxy.predictedCovariance() = *boundParams.covariance();
          trackStateProxy.jacobian() = jacobian;
          trackStateProxy.pathLength() = pathLength;
        } else {
          ACTS_VERBOSE("Detected in-sensitive surface " << surface->geoID());

          // Transport & get curvilinear state instead of bound state
          auto [curvilinearParams, jacobian, pathLength] =
              stepper.curvilinearState(state.stepping, true);

          // Fill the track state
          trackStateProxy.predicted() = curvilinearParams.parameters();
          trackStateProxy.predictedCovariance() =
              *curvilinearParams.covariance();
          trackStateProxy.jacobian() = jacobian;
          trackStateProxy.pathLength() = pathLength;
        }

        // Update state and stepper with material effects
        materialInteractor(surface, state, stepper, fullUpdate);

        // Set the filtered parameter to be the same with predicted parameter
        // @Todo: shall we update the filterd parameter with material effects?
        // But it seems that the smoothing does not like this
        trackStateProxy.filtered() = trackStateProxy.predicted();
        trackStateProxy.filteredCovariance() =
            trackStateProxy.predictedCovariance();

        // We count the processed state
        ++result.processedStates;
      }
      return Result<void>::success();
    }

    /// @brief Track finder actor operation : material interaction
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    /// @tparam stepper_t Type of the stepper
    ///
    /// @param surface The surface where the material interaction happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param updateStage The materal update stage
    ///
    template <typename propagator_state_t, typename stepper_t>
    void materialInteractor(const Surface* surface, propagator_state_t& state,
                            stepper_t& stepper,
                            const MaterialUpdateStage& updateStage) const {
      // Prepare relevant input particle properties
      detail::PointwiseMaterialInteraction interaction(surface, state, stepper);

      // Evaluate the material properties
      if (interaction.evaluateMaterialProperties(state, updateStage)) {
        // Evaluate the material effects
        interaction.evaluatePointwiseMaterialInteraction(multipleScattering,
                                                         energyLoss);

        ACTS_VERBOSE("Material effects on surface: "
                     << surface->geoID() << " at update stage: " << updateStage
                     << " are :");
        ACTS_VERBOSE("eLoss = "
                     << interaction.Eloss << ", "
                     << "variancePhi = " << interaction.variancePhi << ", "
                     << "varianceTheta = " << interaction.varianceTheta << ", "
                     << "varianceQoverP = " << interaction.varianceQoverP);

        // Update the state and stepper with material effects
        interaction.updateState(state, stepper);
      } else {
        ACTS_VERBOSE("No material effects on surface: " << surface->geoID()
                                                        << " at update stage: "
                                                        << updateStage);
      }
    }

    /// Pointer to a logger that is owned by the parent, TrackFinder
    const Logger* m_logger;

    /// Getter for the logger, to support logging macros
    const Logger& logger() const { return *m_logger; }

    /// The track finder updater
    updater_t m_updater;

    /// The track finder smoother
    smoother_t m_smoother;

    /// The Measuremetn calibrator
    calibrator_t m_calibrator;

    /// The Surface beeing
    detail::SurfaceReached targetReached;
  };

  template <typename source_link_t, typename parameters_t>
  class Aborter {
   public:
    /// Broadcast the result_type
    using action_type = Actor<source_link_t, parameters_t>;

    template <typename propagator_state_t, typename stepper_t,
              typename result_t>
    bool operator()(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const result_t& result) const {
      if (!result.result.ok()) {
        return true;
      }
      return false;
    }
  };

 public:
  /// Fit implementation of the foward filter, calls the
  /// the forward filter and backward smoother
  ///
  /// @tparam source_link_t Source link type identifying uncalibrated input
  /// measurements.
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param tfOptions TrackFinderOptions steering the fit
  /// @note The input measurements are given in the form of @c SourceLinks. It's
  /// @c calibrator_t's job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_t, typename start_parameters_t,
            typename parameters_t = BoundParameters,
            typename result_t = Result<TrackFinderResult<source_link_t>>>
  auto findTracks(const std::vector<source_link_t>& sourcelinks,
                  const start_parameters_t& sParameters,
                  const TrackFinderOptions& tfOptions) const
      -> std::enable_if_t<!isDirectNavigator, result_t> {
    static_assert(SourceLinkConcept<source_link_t>,
                  "Source link does not fulfill SourceLinkConcept");

    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_VERBOSE("Preparing " << sourcelinks.size() << " input measurements");
    std::multimap<const Surface*, source_link_t> inputMeasurements;
    for (const auto& sl : sourcelinks) {
      const Surface* srf = &sl.referenceSurface();
      inputMeasurements.emplace(srf, sl);
    }

    // Create the ActionList and AbortList
    using TrackFinderAborter = Aborter<source_link_t, parameters_t>;
    using TrackFinderActor = Actor<source_link_t, parameters_t>;
    using TrackFinderResult = typename TrackFinderActor::result_type;
    using Actors = ActionList<TrackFinderActor>;
    using Aborters = AbortList<TrackFinderAborter>;

    // Create relevant options for the propagation options
    PropagatorOptions<Actors, Aborters> propOptions(tfOptions.geoContext,
                                                    tfOptions.magFieldContext);

    // Catch the actor and set the measurements
    auto& trackFinderActor =
        propOptions.actionList.template get<TrackFinderActor>();
    trackFinderActor.m_logger = m_logger.get();
    trackFinderActor.inputMeasurements = std::move(inputMeasurements);
    trackFinderActor.targetSurface = tfOptions.referenceSurface;
    trackFinderActor.multipleScattering = tfOptions.multipleScattering;
    trackFinderActor.energyLoss = tfOptions.energyLoss;

    // also set logger on updater and smoother
    trackFinderActor.m_updater.m_logger = m_logger;
    trackFinderActor.m_smoother.m_logger = m_logger;

    // Run the fitter
    auto result = m_propagator.template propagate(sParameters, propOptions);

    if (!result.ok()) {
      return result.error();
    }

    const auto& propRes = *result;

    /// Get the result of the fit
    auto trackFinderResult = propRes.template get<TrackFinderResult>();

    /// It could happen that the fit ends in zero processed states.
    /// The result gets meaningless so such case is regarded as fit failure.
    if (trackFinderResult.result.ok() and
        not trackFinderResult.measurementStates) {
      trackFinderResult.result =
          Result<void>(KalmanFitterError::PropagationInVain);
    }

    if (!trackFinderResult.result.ok()) {
      return trackFinderResult.result.error();
    }

    // Return the converted Track
    return m_outputConverter(std::move(trackFinderResult));
  }

  /// Fit implementation of the foward filter, calls the
  /// the forward filter and backward smoother
  ///
  /// @tparam source_link_t Source link type identifying uncalibrated input
  /// measurements.
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  ///
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param tfOptions TrackFinderOptions steering the fit
  /// @param sSequence surface sequence used to initialize a DirectNavigator
  /// @note The input measurements are given in the form of @c SourceLinks. It's
  /// @c calibrator_t's job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_t, typename start_parameters_t,
            typename parameters_t = BoundParameters,
            typename result_t = Result<TrackFinderResult<source_link_t>>>
  auto findTracks(const std::vector<source_link_t>& sourcelinks,
                  const start_parameters_t& sParameters,
                  const TrackFinderOptions& tfOptions,
                  const std::vector<const Surface*>& sSequence) const
      -> std::enable_if_t<isDirectNavigator, result_t> {
    static_assert(SourceLinkConcept<source_link_t>,
                  "Source link does not fulfill SourceLinkConcept");

    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_VERBOSE("Preparing " << sourcelinks.size() << " input measurements");
    std::multimap<const Surface*, source_link_t> inputMeasurements;
    for (const auto& sl : sourcelinks) {
      const Surface* srf = &sl.referenceSurface();
      inputMeasurements.emplace(srf, sl);
    }

    // Create the ActionList and AbortList
    using TrackFinderAborter = Aborter<source_link_t, parameters_t>;
    using TrackFinderActor = Actor<source_link_t, parameters_t>;
    using TrackFinderResult = typename TrackFinderActor::result_type;
    using Actors = ActionList<DirectNavigator::Initializer, TrackFinderActor>;
    using Aborters = AbortList<TrackFinderAborter>;

    // Create relevant options for the propagation options
    PropagatorOptions<Actors, Aborters> propOptions(tfOptions.geoContext,
                                                    tfOptions.magFieldContext);

    // Catch the actor and set the measurements
    auto& trackFinderActor =
        propOptions.actionList.template get<TrackFinderActor>();
    trackFinderActor.m_logger = m_logger.get();
    trackFinderActor.inputMeasurements = std::move(inputMeasurements);
    trackFinderActor.targetSurface = tfOptions.referenceSurface;
    trackFinderActor.multipleScattering = tfOptions.multipleScattering;
    trackFinderActor.energyLoss = tfOptions.energyLoss;

    // also set logger on updater and smoother
    trackFinderActor.m_updater.m_logger = m_logger;
    trackFinderActor.m_smoother.m_logger = m_logger;

    // Set the surface sequence
    auto& dInitializer =
        propOptions.actionList.template get<DirectNavigator::Initializer>();
    dInitializer.surfaceSequence = sSequence;

    // Run the fitter
    auto result = m_propagator.template propagate(sParameters, propOptions);

    if (!result.ok()) {
      return result.error();
    }

    const auto& propRes = *result;

    /// Get the result of the fit
    auto trackFinderResult = propRes.template get<TrackFinderResult>();

    /// It could happen that the fit ends in zero processed states.
    /// The result gets meaningless so such case is regarded as fit failure.
    if (trackFinderResult.result.ok() and
        not trackFinderResult.measurementStates) {
      trackFinderResult.result =
          Result<void>(KalmanFitterError::PropagationInVain);
    }

    if (!trackFinderResult.result.ok()) {
      return trackFinderResult.result.error();
    }

    // Return the converted Track
    return m_outputConverter(std::move(trackFinderResult));
  }
};

}  // namespace Acts

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackStateType.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/PropagatorOptions.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/TrackFitting/detail/KalmanUpdateHelpers.hpp"
#include "Acts/TrackFitting/detail/VoidFitterComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <type_traits>

namespace Acts {

/// Extension struct which holds delegates to customise the KF behavior
template <typename traj_t>
struct KalmanFitterExtensions {
  /// Type alias for track state proxy from trajectory
  using TrackStateProxy = typename traj_t::TrackStateProxy;
  /// Type alias for const track state proxy from trajectory
  using ConstTrackStateProxy = typename traj_t::ConstTrackStateProxy;
  /// Type alias for track parameters from track state proxy
  using Parameters = typename TrackStateProxy::Parameters;

  /// Type alias for measurement calibrator delegate
  using Calibrator =
      Delegate<void(const GeometryContext&, const CalibrationContext&,
                    const SourceLink&, TrackStateProxy)>;

  /// Type alias for Kalman filter update delegate
  using Updater = Delegate<Result<void>(const GeometryContext&, TrackStateProxy,
                                        const Logger&)>;

  /// Type alias for outlier detection delegate
  using OutlierFinder = Delegate<bool(ConstTrackStateProxy)>;

  /// The Calibrator is a dedicated calibration algorithm that allows
  /// to calibrate measurements using track information, this could be
  /// e.g. sagging for wires, module deformations, etc.
  Calibrator calibrator;

  /// The updater incorporates measurement information into the track parameters
  Updater updater;

  /// Determines whether a measurement is supposed to be considered as an
  /// outlier
  OutlierFinder outlierFinder;

  /// Retrieves the associated surface from a source link
  SourceLinkSurfaceAccessor surfaceAccessor;

  /// Default constructor which connects the default void components
  KalmanFitterExtensions() {
    calibrator.template connect<&detail::voidFitterCalibrator<traj_t>>();
    updater.template connect<&detail::voidFitterUpdater<traj_t>>();
    outlierFinder.template connect<&detail::voidOutlierFinder<traj_t>>();
    surfaceAccessor.connect<&detail::voidSurfaceAccessor>();
  }
};

/// Combined options for the Kalman fitter.
///
/// @tparam traj_t The trajectory type
template <typename traj_t>
struct KalmanFitterOptions {
  /// PropagatorOptions with context.
  ///
  /// @param gctx The geometry context for this fit
  /// @param mctx The magnetic context for this fit
  /// @param cctx The calibration context for this fit
  /// @param extensions_ The KF extensions
  /// @param pOptions The plain propagator options
  /// @param tSurface The target surface for the fit
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  /// @param freeToBoundCorrection_ Correction for non-linearity effect during transform from free to bound
  KalmanFitterOptions(const GeometryContext& gctx,
                      const MagneticFieldContext& mctx,
                      std::reference_wrapper<const CalibrationContext> cctx,
                      KalmanFitterExtensions<traj_t> extensions_,
                      const PropagatorPlainOptions& pOptions,
                      const Surface* tSurface = nullptr,
                      bool mScattering = true, bool eLoss = true,
                      const FreeToBoundCorrection& freeToBoundCorrection_ =
                          FreeToBoundCorrection(false))
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        extensions(extensions_),
        propagatorPlainOptions(pOptions),
        targetSurface(tSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        freeToBoundCorrection(freeToBoundCorrection_) {}
  /// Contexts are required and the options must not be default-constructible.
  KalmanFitterOptions() = delete;

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  /// Extensions for calibration and outlier finding
  KalmanFitterExtensions<traj_t> extensions;

  /// The trivial propagator options
  PropagatorPlainOptions propagatorPlainOptions;

  /// The target Surface
  const Surface* targetSurface = nullptr;

  /// Whether to consider multiple scattering
  bool multipleScattering = true;

  /// Whether to consider energy loss
  bool energyLoss = true;

  /// Whether to include non-linear correction during global to local
  /// transformation
  FreeToBoundCorrection freeToBoundCorrection;
};

template <typename traj_t>
struct KalmanFitterResult {
  /// Fitted states that the actor has handled.
  traj_t* fittedStates{nullptr};

  /// This is the index of the 'tip' of the track stored in multitrajectory.
  /// This corresponds to the last measurement state in the multitrajectory.
  /// Since this KF only stores one trajectory, it is unambiguous.
  /// MultiTrajectoryTraits::kInvalid is the start of a trajectory.
  std::size_t lastMeasurementIndex = MultiTrajectoryTraits::kInvalid;

  /// This is the index of the 'tip' of the states stored in multitrajectory.
  /// This corresponds to the last state in the multitrajectory.
  /// Since this KF only stores one trajectory, it is unambiguous.
  /// MultiTrajectoryTraits::kInvalid is the start of a trajectory.
  std::size_t lastTrackIndex = MultiTrajectoryTraits::kInvalid;

  /// The optional Parameters at the provided surface
  std::optional<BoundTrackParameters> fittedParameters;

  /// Counter for states with non-outlier measurements
  std::size_t measurementStates = 0;

  /// Counter for measurements holes
  /// A hole correspond to a surface with an associated detector element with no
  /// associated measurement. Holes are only taken into account if they are
  /// between the first and last measurements.
  std::size_t measurementHoles = 0;

  /// Counter for handled states
  std::size_t processedStates = 0;

  /// Indicator if track fitting has been done
  bool finished = false;

  /// Measurement surfaces without hits
  std::vector<const Surface*> missedActiveSurfaces;

  /// Last encountered error
  Result<void> result{Result<void>::success()};

  /// Path limit aborter
  PathLimitReached pathLimitReached;
};

/// Kalman fitter implementation.
///
/// @tparam propagator_t Type of the propagation class
///
/// The Kalman filter contains an Actor and a Sequencer sub-class.
/// The Sequencer has to be part of the Navigator of the Propagator
/// in order to initialize and provide the measurement surfaces.
///
/// The Actor is part of the Propagation call and does the Kalman update
/// and eventually the smoothing.  Updater, Smoother and Calibrator are
/// given to the Actor for further use:
/// - The Updater is the implemented kalman updater formalism, it
///   runs via a visitor pattern through the measurements.
/// - The Smoother is called at the end of the filtering by the Actor.
///
/// Measurements are not required to be ordered for the KalmanFilter,
/// measurement ordering needs to be figured out by the navigation of
/// the propagator.
///
/// The void components are provided mainly for unit testing.
template <typename propagator_t, typename traj_t>
class KalmanFitter {
  /// The navigator type
  using KalmanNavigator = typename propagator_t::Navigator;

  /// The navigator has DirectNavigator type or not
  static constexpr bool isDirectNavigator =
      std::is_same_v<KalmanNavigator, DirectNavigator>;

 public:
  /// Constructor with propagator and logger
  /// @param pPropagator Propagator instance for track propagation
  /// @param _logger Logger for diagnostic output
  explicit KalmanFitter(propagator_t pPropagator,
                        std::unique_ptr<const Logger> _logger =
                            getDefaultLogger("KalmanFitter", Logging::INFO))
      : m_propagator(std::move(pPropagator)),
        m_logger{std::move(_logger)},
        m_actorLogger{m_logger->cloneWithSuffix("Actor")} {}

 private:
  /// The propagator for the transport and material update
  propagator_t m_propagator;

  /// The logger instance
  std::unique_ptr<const Logger> m_logger;
  std::unique_ptr<const Logger> m_actorLogger;

  const Logger& logger() const { return *m_logger; }

  /// @brief Propagator Actor plugin for the KalmanFilter
  ///
  /// @tparam parameters_t The type of parameters used for "local" parameters.
  /// @tparam calibrator_t The type of calibrator
  /// @tparam outlier_finder_t Type of the outlier finder class
  ///
  /// The KalmanActor does not rely on the measurements to be
  /// sorted along the track.
  template <typename parameters_t>
  class Actor {
   public:
    /// Broadcast the result_type
    using result_type = KalmanFitterResult<traj_t>;

    /// The target surface aboter
    SurfaceReached targetReached{std::numeric_limits<double>::lowest()};

    /// Allows retrieving measurements for a surface
    const std::map<GeometryIdentifier, SourceLink>* inputMeasurements = nullptr;

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// Whether to consider energy loss.
    bool energyLoss = true;

    /// Whether to include non-linear correction during global to local
    /// transformation
    FreeToBoundCorrection freeToBoundCorrection;

    /// Input MultiTrajectory
    std::shared_ptr<traj_t> outputStates;

    KalmanFitterExtensions<traj_t> extensions;

    /// Calibration context for the fit
    const CalibrationContext* calibrationContext{nullptr};

    /// End of world aborter
    EndOfWorldReached endOfWorldReached;

    /// Volume constraint aborter
    VolumeConstraintAborter volumeConstraintAborter;

    /// The logger instance
    const Logger* actorLogger{nullptr};

    /// Logger helper
    const Logger& logger() const { return *actorLogger; }

    /// @brief Kalman actor operation
    ///
    /// @tparam propagator_state_t is the type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void act(propagator_state_t& state, const stepper_t& stepper,
             const navigator_t& navigator, result_type& result,
             const Logger& /*logger*/) const {
      assert(result.fittedStates && "No MultiTrajectory set");

      if (result.finished) {
        return;
      }

      ACTS_VERBOSE("KalmanFitter step at pos: "
                   << stepper.position(state.stepping).transpose()
                   << " dir: " << stepper.direction(state.stepping).transpose()
                   << " momentum: "
                   << stepper.absoluteMomentum(state.stepping));

      // Initialize path limit reached aborter
      if (result.pathLimitReached.internalLimit ==
          std::numeric_limits<double>::max()) {
        detail::setupLoopProtection(state, stepper, result.pathLimitReached,
                                    true, logger());
      }

      // Update:
      // - Waiting for a current surface
      const Surface* surface = navigator.currentSurface(state.navigation);
      if (surface != nullptr) {
        // Check if the surface is in the measurement map
        // -> Get the measurement / calibrate
        // -> Create the predicted state
        // -> Check outlier behavior, if non-outlier:
        // -> Perform the kalman update
        // -> Fill track state information & update stepper information

        ACTS_VERBOSE("Perform " << state.options.direction << " filter step");
        auto res = filter(surface, state, stepper, navigator, result);
        if (!res.ok()) {
          ACTS_ERROR("Error in " << state.options.direction
                                 << " filter: " << res.error());
          result.result = res.error();
        }
      }

      // Finalization:
      // when all track states have been handled or an aborter is triggered
      const bool isTrackComplete =
          result.measurementStates == inputMeasurements->size();
      const bool isEndOfWorldReached =
          endOfWorldReached.checkAbort(state, stepper, navigator, logger());
      const bool isVolumeConstraintReached = volumeConstraintAborter.checkAbort(
          state, stepper, navigator, logger());
      const bool isPathLimitReached = result.pathLimitReached.checkAbort(
          state, stepper, navigator, logger());
      const bool isTargetReached =
          targetReached.checkAbort(state, stepper, navigator, logger());
      if (isTrackComplete || isEndOfWorldReached || isVolumeConstraintReached ||
          isPathLimitReached || isTargetReached) {
        result.finished = true;
      }
    }

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    bool checkAbort(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, const result_type& result,
                    const Logger& /*logger*/) const {
      return !result.result.ok() || result.finished;
    }

    /// @brief Kalman actor operation: update
    ///
    /// @tparam propagator_state_t is the type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param surface The surface where the update happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param result The mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    Result<void> filter(const Surface* surface, propagator_state_t& state,
                        const stepper_t& stepper, const navigator_t& navigator,
                        result_type& result) const {
      const bool precedingMeasurementExists = result.measurementStates > 0;
      const bool surfaceIsSensitive =
          surface->associatedDetectorElement() != nullptr;
      const bool surfaceHasMaterial = surface->surfaceMaterial() != nullptr;

      // Try to find the surface in the measurement surfaces
      const auto sourceLinkIt = inputMeasurements->find(surface->geometryId());
      if (sourceLinkIt != inputMeasurements->end()) {
        // Screen output message
        ACTS_VERBOSE("Measurement surface " << surface->geometryId()
                                            << " detected.");
        // Transport the covariance to the surface
        stepper.transportCovarianceToBound(state.stepping, *surface,
                                           freeToBoundCorrection);

        // Update state and stepper with pre material effects
        materialInteractor(surface, state, stepper, navigator,
                           MaterialUpdateStage::PreUpdate);

        // do the kalman update (no need to perform covTransport here, hence no
        // point in performing globalToLocal correction)
        const auto trackStateProxyRes = detail::kalmanHandleMeasurement(
            *calibrationContext, state, stepper, extensions, *surface,
            sourceLinkIt->second, *result.fittedStates, result.lastTrackIndex,
            false, logger());

        if (!trackStateProxyRes.ok()) {
          return trackStateProxyRes.error();
        }

        const auto& trackStateProxy = *trackStateProxyRes;
        result.lastTrackIndex = trackStateProxy.index();

        // Update the stepper if it is not an outlier
        if (trackStateProxy.typeFlags().test(TrackStateFlag::MeasurementFlag)) {
          // Update the stepping state with filtered parameters
          ACTS_VERBOSE("Filtering step successful, updated parameters are:\n"
                       << trackStateProxy.filtered().transpose());
          // update stepping state using filtered parameters after kalman
          stepper.update(state.stepping,
                         MultiTrajectoryHelpers::freeFiltered(
                             state.options.geoContext, trackStateProxy),
                         trackStateProxy.filtered(),
                         trackStateProxy.filteredCovariance(), *surface);
          // We count the state with measurement
          ++result.measurementStates;
        }

        // Update state and stepper with post material effects
        materialInteractor(surface, state, stepper, navigator,
                           MaterialUpdateStage::PostUpdate);
        // We count the processed state
        ++result.processedStates;
        // Update the number of holes count only when encountering a
        // measurement
        result.measurementHoles = result.missedActiveSurfaces.size();
        // Since we encountered a measurement update the lastMeasurementIndex to
        // the lastTrackIndex.
        result.lastMeasurementIndex = result.lastTrackIndex;

      } else if ((precedingMeasurementExists && surfaceIsSensitive) ||
                 surfaceHasMaterial) {
        // We only create track states here if there is already measurement
        // detected or if the surface has material (no holes before the first
        // measurement)
        const auto trackStateProxyRes = detail::kalmanHandleNoMeasurement(
            state.stepping, stepper, *surface, *result.fittedStates,
            result.lastTrackIndex, true, logger(), precedingMeasurementExists,
            freeToBoundCorrection);

        if (!trackStateProxyRes.ok()) {
          return trackStateProxyRes.error();
        }

        const auto& trackStateProxy = *trackStateProxyRes;
        result.lastTrackIndex = trackStateProxy.index();

        if (trackStateProxy.typeFlags().test(TrackStateFlag::HoleFlag)) {
          // Count the missed surface
          result.missedActiveSurfaces.push_back(surface);
        }

        ++result.processedStates;

        // Update state and stepper with (possible) material effects
        materialInteractor(surface, state, stepper, navigator,
                           MaterialUpdateStage::FullUpdate);
      }

      return Result<void>::success();
    }

    /// @brief Kalman actor operation: material interaction
    ///
    /// @tparam propagator_state_t is the type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param surface The surface where the material interaction happens
    /// @param state The mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param updateStage The material update stage
    ///
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    void materialInteractor(const Surface* surface, propagator_state_t& state,
                            const stepper_t& stepper,
                            const navigator_t& navigator,
                            const MaterialUpdateStage& updateStage) const {
      // Protect against null surface
      if (surface == nullptr) {
        ACTS_VERBOSE(
            "Surface is nullptr. Cannot be used for material interaction");
        return;
      }

      if (surface->surfaceMaterial() == nullptr) {
        ACTS_VERBOSE("No material on surface: " << surface->geometryId());
        return;
      }

      // Prepare relevant input particle properties
      detail::PointwiseMaterialInteraction interaction(surface, state, stepper);

      if (!interaction.evaluateMaterialSlab(state, navigator, updateStage)) {
        ACTS_VERBOSE("No material on surface after evaluation: "
                     << surface->geometryId());
        return;
      }

      // Evaluate the material effects
      interaction.evaluatePointwiseMaterialInteraction(multipleScattering,
                                                       energyLoss);

      // Screen out material effects info
      ACTS_VERBOSE("Material effects on surface: " << surface->geometryId()
                                                   << " at update stage: "
                                                   << updateStage << " are :");
      ACTS_VERBOSE("eLoss = "
                   << interaction.Eloss << ", "
                   << "variancePhi = " << interaction.variancePhi << ", "
                   << "varianceTheta = " << interaction.varianceTheta << ", "
                   << "varianceQoverP = " << interaction.varianceQoverP);

      // Update the state and stepper with material effects
      interaction.updateState(state, stepper, addNoise);
    }
  };

 public:
  /// Fit implementation of the forward filter, calls the
  /// the filter and smoother/reversed filter
  ///
  /// @tparam source_link_iterator_t Iterator type used to pass source links
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  /// @tparam track_container_t Type of the track container
  ///
  /// @param it Begin iterator for the fittable uncalibrated measurements
  /// @param end End iterator for the fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param kfOptions KalmanOptions steering the fit
  /// @param trackContainer Input track container storage to append into
  /// @note The input measurements are given in the form of @c SourceLink s.
  /// It's the calibrators job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_iterator_t, typename start_parameters_t,
            typename parameters_t = BoundTrackParameters,
            TrackContainerFrontend track_container_t>
  Result<typename track_container_t::TrackProxy> fit(
      source_link_iterator_t it, source_link_iterator_t end,
      const start_parameters_t& sParameters,
      const KalmanFitterOptions<traj_t>& kfOptions,
      track_container_t& trackContainer) const
    requires(!isDirectNavigator)
  {
    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyway, so the map can own them.
    ACTS_VERBOSE("Preparing " << std::distance(it, end)
                              << " input measurements");
    std::map<GeometryIdentifier, SourceLink> inputMeasurements;
    // for (const auto& sl : sourceLinks) {
    for (; it != end; ++it) {
      SourceLink sl = *it;
      const Surface* surface = kfOptions.extensions.surfaceAccessor(sl);
      // @TODO: This can probably change over to surface pointers as keys
      auto geoId = surface->geometryId();
      inputMeasurements.emplace(geoId, std::move(sl));
    }

    // Create the ActorList
    using KalmanActor = Actor<parameters_t>;

    using KalmanResult = typename KalmanActor::result_type;
    using Actors = ActorList<KalmanActor>;
    using PropagatorOptions = typename propagator_t::template Options<Actors>;

    // Create relevant options for the propagation options
    PropagatorOptions propagatorOptions(kfOptions.geoContext,
                                        kfOptions.magFieldContext);

    // Set the trivial propagator options
    propagatorOptions.setPlainOptions(kfOptions.propagatorPlainOptions);

    // Add the measurement surface as external surface to navigator.
    // We will try to hit those surface by ignoring boundary checks.
    for (const auto& [surfaceId, _] : inputMeasurements) {
      propagatorOptions.navigation.insertExternalSurface(surfaceId);
    }

    // Catch the actor and set the measurements
    auto& kalmanActor = propagatorOptions.actorList.template get<KalmanActor>();
    kalmanActor.inputMeasurements = &inputMeasurements;
    kalmanActor.targetReached.surface = kfOptions.targetSurface;
    kalmanActor.multipleScattering = kfOptions.multipleScattering;
    kalmanActor.energyLoss = kfOptions.energyLoss;
    kalmanActor.freeToBoundCorrection = kfOptions.freeToBoundCorrection;
    kalmanActor.calibrationContext = &kfOptions.calibrationContext.get();
    kalmanActor.extensions = kfOptions.extensions;
    kalmanActor.actorLogger = m_actorLogger.get();

    return fit_impl<start_parameters_t, PropagatorOptions, KalmanResult,
                    track_container_t>(sParameters, propagatorOptions,
                                       trackContainer);
  }

  /// Fit implementation of the forward filter, calls the
  /// the filter and smoother/reversed filter
  ///
  /// @tparam source_link_iterator_t Iterator type used to pass source links
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam parameters_t Type of parameters used for local parameters
  /// @tparam track_container_t Type of the track container
  ///
  /// @param it Begin iterator for the fittable uncalibrated measurements
  /// @param end End iterator for the fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param kfOptions KalmanOptions steering the fit
  /// @param sSequence surface sequence used to initialize a DirectNavigator
  /// @param trackContainer Input track container storage to append into
  /// @note The input measurements are given in the form of @c SourceLinks.
  /// It's
  /// @c calibrator_t's job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_iterator_t, typename start_parameters_t,
            typename parameters_t = BoundTrackParameters,
            TrackContainerFrontend track_container_t>
  Result<typename track_container_t::TrackProxy> fit(
      source_link_iterator_t it, source_link_iterator_t end,
      const start_parameters_t& sParameters,
      const KalmanFitterOptions<traj_t>& kfOptions,
      const std::vector<const Surface*>& sSequence,
      track_container_t& trackContainer) const
    requires(isDirectNavigator)
  {
    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyway, so the map can own them.
    ACTS_VERBOSE("Preparing " << std::distance(it, end)
                              << " input measurements");
    std::map<GeometryIdentifier, SourceLink> inputMeasurements;
    for (; it != end; ++it) {
      SourceLink sl = *it;
      const Surface* surface = kfOptions.extensions.surfaceAccessor(sl);
      // @TODO: This can probably change over to surface pointers as keys
      auto geoId = surface->geometryId();
      inputMeasurements.emplace(geoId, std::move(sl));
    }

    // Create the ActorList
    using KalmanActor = Actor<parameters_t>;

    using KalmanResult = typename KalmanActor::result_type;
    using Actors = ActorList<KalmanActor>;
    using PropagatorOptions = typename propagator_t::template Options<Actors>;

    // Create relevant options for the propagation options
    PropagatorOptions propagatorOptions(kfOptions.geoContext,
                                        kfOptions.magFieldContext);

    // Set the trivial propagator options
    propagatorOptions.setPlainOptions(kfOptions.propagatorPlainOptions);

    // Catch the actor and set the measurements
    auto& kalmanActor = propagatorOptions.actorList.template get<KalmanActor>();
    kalmanActor.inputMeasurements = &inputMeasurements;
    kalmanActor.targetReached.surface = kfOptions.targetSurface;
    kalmanActor.multipleScattering = kfOptions.multipleScattering;
    kalmanActor.energyLoss = kfOptions.energyLoss;
    kalmanActor.freeToBoundCorrection = kfOptions.freeToBoundCorrection;
    kalmanActor.calibrationContext = &kfOptions.calibrationContext.get();
    kalmanActor.extensions = kfOptions.extensions;
    kalmanActor.actorLogger = m_actorLogger.get();

    // Set the surface sequence
    propagatorOptions.navigation.surfaces = sSequence;

    return fit_impl<start_parameters_t, PropagatorOptions, KalmanResult,
                    track_container_t>(sParameters, propagatorOptions,
                                       trackContainer);
  }

 private:
  /// Common fit implementation
  ///
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam actor_list_t Type of the actor list
  /// @tparam aborter_list_t Type of the abort list
  /// @tparam kalman_result_t Type of the KF result
  /// @tparam track_container_t Type of the track container
  ///
  /// @param sParameters The initial track parameters
  /// @param propagatorOptions The Propagator Options
  /// @param trackContainer Input track container storage to append into
  ///
  /// @return the output as an output track
  template <typename start_parameters_t, typename propagator_options_t,
            typename kalman_result_t, TrackContainerFrontend track_container_t>
  auto fit_impl(const start_parameters_t& sParameters,
                const propagator_options_t& propagatorOptions,
                track_container_t& trackContainer) const
      -> Result<typename track_container_t::TrackProxy> {
    auto propagatorState = m_propagator.makeState(propagatorOptions);

    auto propagatorInitResult =
        m_propagator.initialize(propagatorState, sParameters);
    if (!propagatorInitResult.ok()) {
      ACTS_ERROR("Propagation initialization failed: "
                 << propagatorInitResult.error());
      return propagatorInitResult.error();
    }

    auto& kalmanResult =
        propagatorState.template get<KalmanFitterResult<traj_t>>();
    kalmanResult.fittedStates = &trackContainer.trackStateContainer();

    // Run the fitter
    auto result = m_propagator.propagate(propagatorState);

    if (!result.ok()) {
      ACTS_ERROR("Propagation failed: " << result.error());
      return result.error();
    }

    // It could happen that the fit ends in zero measurement states.
    // The result gets meaningless so such case is regarded as fit failure.
    if (kalmanResult.result.ok() && !kalmanResult.measurementStates) {
      kalmanResult.result = Result<void>(KalmanFitterError::NoMeasurementFound);
    }

    if (!kalmanResult.result.ok()) {
      ACTS_ERROR("KalmanFilter failed: "
                 << kalmanResult.result.error() << ", "
                 << kalmanResult.result.error().message());
      return kalmanResult.result.error();
    }

    auto track = trackContainer.makeTrack();
    track.tipIndex() = kalmanResult.lastMeasurementIndex;
    if (kalmanResult.fittedParameters) {
      const auto& params = kalmanResult.fittedParameters.value();
      track.parameters() = params.parameters();
      track.covariance() = params.covariance().value();
      track.setReferenceSurface(params.referenceSurface().getSharedPtr());
    }

    calculateTrackQuantities(track);

    return track;
  }
};

}  // namespace Acts

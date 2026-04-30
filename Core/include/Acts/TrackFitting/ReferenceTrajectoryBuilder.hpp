// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackContainerFrontendConcept.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/PropagatorOptions.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <vector>

namespace Acts::Experimental {

/// @tparam traj_t The trajectory type
template <typename traj_t>
struct ReferenceTrajectoryBuilderOptions {
  /// PropagatorOptions with context.
  ///
  /// @param gctx The geometry context
  /// @param mctx The magnetic context
  /// @param pOptions The plain propagator options
  /// @param tSurface The target surface for the fit
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  /// @param freeToBoundCorrection_ Correction for non-linearity effect during transform from free to bound
  ReferenceTrajectoryBuilderOptions(
      const GeometryContext& gctx, const MagneticFieldContext& mctx,
      const PropagatorPlainOptions& pOptions, const Surface* tSurface = nullptr,
      bool mScattering = true, bool eLoss = true,
      const FreeToBoundCorrection& freeToBoundCorrection_ =
          FreeToBoundCorrection(false))
      : geoContext(gctx),
        magFieldContext(mctx),
        propagatorPlainOptions(pOptions),
        referenceSurface(tSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        freeToBoundCorrection(freeToBoundCorrection_) {}

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;

  /// The trivial propagator options
  PropagatorPlainOptions propagatorPlainOptions;

  /// The reference surface
  const Surface* referenceSurface = nullptr;

  /// Whether to consider multiple scattering
  bool multipleScattering = true;

  /// Whether to consider energy loss
  bool energyLoss = true;

  /// Whether to include non-linear correction during global to local
  /// transformation
  FreeToBoundCorrection freeToBoundCorrection;
};

template <typename traj_t>
struct ReferenceTrajectoryBuilderResult {
  traj_t* trajectory{nullptr};

  std::size_t lastTrackStateIndex = kTrackIndexInvalid;

  std::optional<BoundTrackParameters> referenceParameters;

  bool finished = false;

  /// Path limit aborter
  PathLimitReached pathLimitReached;
};

/// Reference trajectory implementation.
template <typename propagator_t, typename traj_t>
class ReferenceTrajectoryBuilder {
  /// The navigator type
  using NavigatorType = typename propagator_t::Navigator;

  /// The navigator has DirectNavigator type or not
  static constexpr bool isDirectNavigator =
      std::is_same_v<NavigatorType, DirectNavigator>;

  using ResultType = ReferenceTrajectoryBuilderResult<traj_t>;

 public:
  using Options = ReferenceTrajectoryBuilderOptions<traj_t>;

  /// Type alias for track state proxy from trajectory
  using TrackStateProxy = typename traj_t::TrackStateProxy;

  using ConstTrackStateProxy = typename traj_t::ConstTrackStateProxy;

  /// Constructor with propagator and logger
  /// @param pPropagator Propagator instance for track propagation
  /// @param _logger Logger for diagnostic output
  explicit ReferenceTrajectoryBuilder(
      propagator_t pPropagator,
      std::unique_ptr<const Logger> _logger =
          getDefaultLogger("ReferenceTrajectoryBuilder", Logging::INFO))
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

  class Actor {
   public:
    /// Broadcast the result_type
    using result_type = ResultType;

    /// The target surface aborter
    SurfaceReached targetReached{std::numeric_limits<double>::lowest()};

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// Whether to consider energy loss.
    bool energyLoss = true;

    /// Whether to include non-linear correction during global to local
    /// transformation
    FreeToBoundCorrection freeToBoundCorrection;

    /// Input MultiTrajectory
    std::shared_ptr<traj_t> outputStates;

    /// End of world aborter
    EndOfWorldReached endOfWorldReached;

    /// Volume constraint aborter
    VolumeConstraintAborter volumeConstraintAborter;

    /// The logger instance
    const Logger* actorLogger{nullptr};

    /// Logger helper
    const Logger& logger() const { return *actorLogger; }

    /// Actor operation
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
    Result<void> act(propagator_state_t& state, const stepper_t& stepper,
                     const navigator_t& navigator, result_type& result,
                     const Logger& /*logger*/) const {
      assert(result.trajectory && "No MultiTrajectory set");

      if (result.finished) {
        return Result<void>::success();
      }

      ACTS_VERBOSE("ReferenceTrajectory step at pos: "
                   << stepper.position(state.stepping).transpose()
                   << " dir: " << stepper.direction(state.stepping).transpose()
                   << " momentum: "
                   << stepper.absoluteMomentum(state.stepping));

      if (result.pathLimitReached.internalLimit ==
          std::numeric_limits<double>::max()) {
        detail::setupLoopProtection(state, stepper, result.pathLimitReached,
                                    true, logger());
      }

      if (const Surface* surface = navigator.currentSurface(state.navigation);
          surface != nullptr) {
        ACTS_VERBOSE("Handle Surface " << surface->geometryId() << " "
                                       << state.options.direction);
        auto res = handleSurface(*surface, state, stepper, navigator, result);
        if (!res.ok()) {
          ACTS_DEBUG("Error in " << state.options.direction
                                 << " filter: " << res.error());
          return res.error();
        }
      }

      const bool isEndOfWorldReached =
          endOfWorldReached.checkAbort(state, stepper, navigator, logger());
      const bool isVolumeConstraintReached = volumeConstraintAborter.checkAbort(
          state, stepper, navigator, logger());
      const bool isPathLimitReached = result.pathLimitReached.checkAbort(
          state, stepper, navigator, logger());
      const bool isTargetReached =
          targetReached.checkAbort(state, stepper, navigator, logger());
      if (isEndOfWorldReached || isVolumeConstraintReached ||
          isPathLimitReached || isTargetReached) {
        ACTS_VERBOSE(
            "Finalizing reference trajectory: "
            << (isEndOfWorldReached ? "end of world reached; " : "")
            << (isVolumeConstraintReached ? "volume constraint reached; " : "")
            << (isPathLimitReached ? "path limit reached; " : "")
            << (isTargetReached ? "target surface reached; " : ""));

        if (isTargetReached) {
          ACTS_VERBOSE("Setting parameters at target surface");

          auto res = stepper.boundState(state.stepping, *targetReached.surface);
          if (!res.ok()) {
            ACTS_DEBUG("Error while acquiring bound state for target surface: "
                       << res.error() << " " << res.error().message());
            return res.error();
          } else {
            const auto& [boundParams, jacobian, pathLength] = *res;
            result.referenceParameters = boundParams;
          }
        }

        result.finished = true;
      }

      return Result<void>::success();
    }

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    bool checkAbort(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, const result_type& result,
                    const Logger& /*logger*/) const {
      return result.finished;
    }

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    Result<void> handleSurface(const Surface& surface,
                               propagator_state_t& state,
                               const stepper_t& stepper,
                               const navigator_t& navigator,
                               result_type& result) const {
      stepper.transportCovarianceToBound(state.stepping, surface,
                                         freeToBoundCorrection);

      detail::performMaterialInteraction(
          state, stepper, surface,
          detail::determineMaterialUpdateMode(state, navigator,
                                              MaterialUpdateMode::PreUpdate),
          NoiseUpdateMode::addNoise, multipleScattering, energyLoss, logger());

      TrackStatePropMask mask =
          TrackStatePropMask::Predicted | TrackStatePropMask::Jacobian;
      TrackStateProxy trackStateProxy =
          result.trajectory->makeTrackState(mask, result.lastTrackStateIndex);

      ConstTrackStateProxy trackStateProxyConst{trackStateProxy};

      trackStateProxy.setReferenceSurface(surface.getSharedPtr());
      auto res = stepper.boundState(state.stepping, surface, false,
                                    freeToBoundCorrection);
      if (!res.ok()) {
        ACTS_DEBUG("Propagate to surface " << surface.geometryId()
                                           << " failed: " << res.error());
        return res.error();
      }
      const auto& [boundParams, jacobian, pathLength] = *res;

      trackStateProxy.predicted() = boundParams.parameters();
      trackStateProxy.predictedCovariance() = state.stepping.cov;
      trackStateProxy.jacobian() = jacobian;
      trackStateProxy.pathLength() = pathLength;

      auto typeFlags = trackStateProxy.typeFlags();
      typeFlags.setHasParameters();
      if (surface.surfaceMaterial() != nullptr) {
        typeFlags.setHasMaterial();
      }

      result.lastTrackStateIndex = trackStateProxy.index();

      detail::performMaterialInteraction(
          state, stepper, surface,
          detail::determineMaterialUpdateMode(state, navigator,
                                              MaterialUpdateMode::PostUpdate),
          NoiseUpdateMode::addNoise, multipleScattering, energyLoss, logger());

      return Result<void>::success();
    }
  };

 public:
  template <TrackContainerFrontend track_container_t>
  Result<typename track_container_t::TrackProxy> build(
      const BoundTrackParameters& sParameters, const Options& actorOptions,
      track_container_t& trackContainer) const {
    auto propagatorOptions = makePropagatorOptions(
        actorOptions, nullptr, actorOptions.referenceSurface);
    return buildImpl(sParameters, propagatorOptions, trackContainer);
  }

  template <TrackContainerFrontend track_container_t>
  Result<typename track_container_t::TrackProxy> build(
      const BoundTrackParameters& sParameters, const Options& actorOptions,
      const std::vector<const Surface*>& sSequence,
      track_container_t& trackContainer) const
    requires(isDirectNavigator)
  {
    auto propagatorOptions = makePropagatorOptions(
        actorOptions, &sSequence, actorOptions.referenceSurface);
    return buildImpl(sParameters, propagatorOptions, trackContainer);
  }

  template <typename track_proxy_t, typename source_link_range_t>
  void attachSourceLinks(
      track_proxy_t trackProxy, source_link_range_t sourceLinkRange,
      const SourceLinkSurfaceAccessor& surfaceAccessor) const {
    const std::size_t nMeasurements = std::ranges::distance(sourceLinkRange);

    ACTS_VERBOSE("Preparing " << nMeasurements << " input measurements");
    std::unordered_map<const Surface*, SourceLink> inputMeasurements;
    for (auto&& sl : sourceLinkRange) {
      const Surface* surface = surfaceAccessor(sl);
      inputMeasurements.try_emplace(surface, std::move(sl));
    }

    for (auto trackState : trackProxy.trackStates()) {
      if (!trackState.hasReferenceSurface()) {
        continue;
      }
      const Surface& surface = trackState.referenceSurface();

      if (!surface.isSensitive()) {
        continue;
      }

      auto typeFlagsMap = trackState.typeFlags();

      if (const auto it = inputMeasurements.find(&surface);
          it == inputMeasurements.end()) {
        typeFlagsMap.setIsHole();
      } else {
        auto sourceLink = it->second;

        trackState.setUncalibratedSourceLink(std::move(sourceLink));
        typeFlagsMap.setHasMeasurement();
      }
    }
  }

  using Calibrator =
      Delegate<void(const GeometryContext&, const CalibrationContext&,
                    const SourceLink&, TrackStateProxy)>;

  template <typename track_proxy_t>
  void calibrateMeasurements(const GeometryContext& geoContext,
                             const CalibrationContext& calibrationContext,
                             track_proxy_t trackProxy,
                             const Calibrator& calibrator) const {
    for (auto trackState : trackProxy.trackStates()) {
      if (!trackState.typeFlags().hasMeasurement()) {
        continue;
      }

      trackState.addComponents(TrackStatePropMask::Calibrated);

      const SourceLink& sourceLink = trackState.getUncalibratedSourceLink();
      calibrator(geoContext, calibrationContext, sourceLink, trackState);
    }
  }

  using Updater = Delegate<Result<void>(const GeometryContext&, TrackStateProxy,
                                        const Logger&)>;

  template <typename track_proxy_t>
  Result<void> filter(const GeometryContext& geoContext,
                      track_proxy_t trackProxy, const Updater& updater) const {
    BoundVector lastDelta = BoundVector::Zero();

    for (auto trackState : trackProxy.trackStates()) {
      if (!trackState.typeFlags().hasMeasurement()) {
        trackState.shareFrom(trackState, TrackStatePropMask::Predicted,
                             TrackStatePropMask::Filtered);
        continue;
      }

      trackState.addComponents(TrackStatePropMask::Filtered);

      // Transport the last delta to the current surface using the Jacobian of
      // the track state
      trackState.predicted() += trackState.jacobian() * lastDelta;

      auto res = updater(geoContext, trackState, logger());
      if (!res.ok()) {
        ACTS_DEBUG("Error in filter: " << res.error());
        return res.error();
      }

      // Update the last delta with the change in parameters after the update
      lastDelta = trackState.filtered() - trackState.predicted();
    }

    return Result<void>::success();
  }

 private:
  auto makePropagatorOptions(const Options& actorOptions,
                             const std::vector<const Surface*>* sSequence,
                             const Surface* targetSurface) const {
    using Actors = ActorList<Actor>;
    using PropagatorOptions = typename propagator_t::template Options<Actors>;

    PropagatorOptions propagatorOptions(actorOptions.geoContext,
                                        actorOptions.magFieldContext);
    propagatorOptions.setPlainOptions(actorOptions.propagatorPlainOptions);

    if constexpr (!isDirectNavigator) {
      if (sSequence != nullptr) {
        for (const Surface* surface : *sSequence) {
          propagatorOptions.navigation.appendExternalSurface(*surface);
        }
      }
    } else {
      assert(sSequence != nullptr &&
             "DirectNavigator requires a surface sequence for "
             "ReferenceTrajectory");
      propagatorOptions.navigation.externalSurfaces = *sSequence;
    }

    auto& actor = propagatorOptions.actorList.template get<Actor>();
    actor.targetReached.surface = targetSurface;
    actor.multipleScattering = actorOptions.multipleScattering;
    actor.energyLoss = actorOptions.energyLoss;
    actor.freeToBoundCorrection = actorOptions.freeToBoundCorrection;
    actor.actorLogger = m_actorLogger.get();

    return propagatorOptions;
  }

  template <typename propagator_options_t,
            TrackContainerFrontend track_container_t>
  auto buildImpl(const BoundTrackParameters& sParameters,
                 const propagator_options_t& propagatorOptions,
                 track_container_t& trackContainer) const
      -> Result<typename track_container_t::TrackProxy> {
    auto propagatorState = m_propagator.makeState(propagatorOptions);

    auto propagatorInitResult =
        m_propagator.initialize(propagatorState, sParameters);
    if (!propagatorInitResult.ok()) {
      ACTS_DEBUG("Propagation initialization failed: "
                 << propagatorInitResult.error());
      return propagatorInitResult.error();
    }

    auto& actorResult = propagatorState.template get<ResultType>();
    actorResult.trajectory = &trackContainer.trackStateContainer();

    auto result = m_propagator.propagate(propagatorState);

    if (!result.ok()) {
      ACTS_DEBUG("Propagation failed: " << result.error());
      return result.error();
    }

    auto track = trackContainer.makeTrack();
    track.tipIndex() = actorResult.lastTrackStateIndex;
    if (actorResult.referenceParameters.has_value()) {
      const auto& params = *actorResult.referenceParameters;
      track.parameters() = params.parameters();
      track.covariance() = params.covariance().value();
      track.setReferenceSurface(params.referenceSurface().getSharedPtr());
    }

    track.linkForward();

    return track;
  }
};

}  // namespace Acts::Experimental

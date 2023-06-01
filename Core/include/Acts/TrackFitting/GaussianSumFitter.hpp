// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackHelpers.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/MultiStepperAborters.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/GsfActor.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <fstream>

namespace Acts {

namespace detail {

/// Type trait to identify if a type is a MultiComponentBoundTrackParameters and
/// to inspect its charge representation if not TODO this probably gives an ugly
/// error message if detectCharge does not compile
template <typename T>
struct IsMultiComponentBoundParameters : public std::false_type {
  template <template <class> class U, class V>
  static auto detectCharge(const U<V>& /*unused*/) {
    return V{};
  }

  using Charge = decltype(detectCharge(std::declval<T>()));
};

template <typename T>
struct IsMultiComponentBoundParameters<MultiComponentBoundTrackParameters<T>>
    : public std::true_type {};

}  // namespace detail

namespace Experimental {

/// Gaussian Sum Fitter implementation.
/// @tparam propagator_t The propagator type on which the algorithm is built on
/// @tparam bethe_heitler_approx_t The type of the Bethe-Heitler-Approximation
/// @tparam traj_t The MultiTrajectory type (backend)
///
/// @note This GSF implementation tries to be as compatible to the KalmanFitter
/// as possible. However, strict compatibility is not garantueed.
/// @note Currently there is no possibility to export the states of the
/// individual components from the GSF, the only information returned in the
/// MultiTrajectory are the means of the states. Therefore, also NO dedicated
/// component smoothing is performed as described e.g. by R. Fruewirth.
template <typename propagator_t, typename bethe_heitler_approx_t,
          typename traj_t>
struct GaussianSumFitter {
  GaussianSumFitter(propagator_t&& propagator, bethe_heitler_approx_t&& bha,
                    std::unique_ptr<const Logger> _logger =
                        getDefaultLogger("GSF", Logging::INFO))
      : m_propagator(std::move(propagator)),
        m_betheHeitlerApproximation(std::move(bha)),
        m_logger{std::move(_logger)},
        m_actorLogger(m_logger->cloneWithSuffix("Actor")) {}

  /// The propagator instance used by the fit function
  propagator_t m_propagator;

  /// The fitter holds the instance of the bethe heitler approx
  bethe_heitler_approx_t m_betheHeitlerApproximation;

  /// The logger
  std::unique_ptr<const Logger> m_logger;
  std::unique_ptr<const Logger> m_actorLogger;

  const Logger& logger() const { return *m_logger; }

  /// The navigator type
  using GsfNavigator = typename propagator_t::Navigator;

  /// The actor type
  using GsfActor = detail::GsfActor<bethe_heitler_approx_t, traj_t>;

  /// @brief The fit function for the Direct navigator
  template <typename source_link_it_t, typename start_parameters_t,
            typename track_container_t, template <typename> class holder_t>
  auto fit(source_link_it_t begin, source_link_it_t end,
           const start_parameters_t& sParameters,
           const GsfOptions<traj_t>& options,
           const std::vector<const Surface*>& sSequence,
           TrackContainer<track_container_t, traj_t, holder_t>& trackContainer)
      const {
    // Check if we have the correct navigator
    static_assert(
        std::is_same_v<DirectNavigator, typename propagator_t::Navigator>);

    // Initialize the forward propagation with the DirectNavigator
    auto fwdPropInitializer = [&sSequence, this](const auto& opts) {
      using Actors = ActionList<GsfActor, DirectNavigator::Initializer>;
      using Aborters = AbortList<>;

      PropagatorOptions<Actors, Aborters> propOptions(opts.geoContext,
                                                      opts.magFieldContext);

      propOptions.setPlainOptions(opts.propagatorPlainOptions);

      propOptions.actionList.template get<DirectNavigator::Initializer>()
          .navSurfaces = sSequence;
      propOptions.actionList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_betheHeitlerApproximation;

      return propOptions;
    };

    // Initialize the backward propagation with the DirectNavigator
    auto bwdPropInitializer = [&sSequence, this](const auto& opts) {
      using Actors = ActionList<GsfActor, Acts::detail::FinalStateCollector,
                                DirectNavigator::Initializer>;
      using Aborters = AbortList<>;

      std::vector<const Surface*> backwardSequence(
          std::next(sSequence.rbegin()), sSequence.rend());
      backwardSequence.push_back(opts.referenceSurface);

      PropagatorOptions<Actors, Aborters> propOptions(opts.geoContext,
                                                      opts.magFieldContext);

      propOptions.setPlainOptions(opts.propagatorPlainOptions);

      propOptions.actionList.template get<DirectNavigator::Initializer>()
          .navSurfaces = std::move(backwardSequence);
      propOptions.actionList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_betheHeitlerApproximation;

      return propOptions;
    };

    return fit_impl(begin, end, sParameters, options, fwdPropInitializer,
                    bwdPropInitializer, trackContainer);
  }

  /// @brief The fit function for the standard navigator
  template <typename source_link_it_t, typename start_parameters_t,
            typename track_container_t, template <typename> class holder_t>
  auto fit(source_link_it_t begin, source_link_it_t end,
           const start_parameters_t& sParameters,
           const GsfOptions<traj_t>& options,
           TrackContainer<track_container_t, traj_t, holder_t>& trackContainer)
      const {
    // Check if we have the correct navigator
    static_assert(std::is_same_v<Navigator, typename propagator_t::Navigator>);

    // Initialize the forward propagation with the DirectNavigator
    auto fwdPropInitializer = [this](const auto& opts) {
      using Actors = ActionList<GsfActor>;
      using Aborters = AbortList<EndOfWorldReached>;

      PropagatorOptions<Actors, Aborters> propOptions(opts.geoContext,
                                                      opts.magFieldContext);
      propOptions.setPlainOptions(opts.propagatorPlainOptions);
      propOptions.actionList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_betheHeitlerApproximation;

      return propOptions;
    };

    // Initialize the backward propagation with the DirectNavigator
    auto bwdPropInitializer = [this](const auto& opts) {
      using Actors = ActionList<GsfActor, Acts::detail::FinalStateCollector>;
      using Aborters = AbortList<EndOfWorldReached>;

      PropagatorOptions<Actors, Aborters> propOptions(opts.geoContext,
                                                      opts.magFieldContext);

      propOptions.setPlainOptions(opts.propagatorPlainOptions);

      propOptions.actionList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_betheHeitlerApproximation;

      return propOptions;
    };

    return fit_impl(begin, end, sParameters, options, fwdPropInitializer,
                    bwdPropInitializer, trackContainer);
  }

  /// The generic implementation of the fit function.
  /// TODO check what this function does with the referenceSurface is e.g. the
  /// first measuerementSurface
  template <typename source_link_it_t, typename start_parameters_t,
            typename fwd_prop_initializer_t, typename bwd_prop_initializer_t,
            typename track_container_t, template <typename> class holder_t>
  Acts::Result<
      typename TrackContainer<track_container_t, traj_t, holder_t>::TrackProxy>
  fit_impl(source_link_it_t begin, source_link_it_t end,
           const start_parameters_t& sParameters,
           const GsfOptions<traj_t>& options,
           const fwd_prop_initializer_t& fwdPropInitializer,
           const bwd_prop_initializer_t& bwdPropInitializer,
           TrackContainer<track_container_t, traj_t, holder_t>& trackContainer)
      const {
    // return or abort utility
    auto return_error_or_abort = [&](auto error) {
      if (options.abortOnError) {
        std::abort();
      }
      return error;
    };

    // Define directions based on input propagation direction. This way we can
    // refer to 'forward' and 'backward' regardless of the actual direction.
    const auto gsfForward = options.propagatorPlainOptions.direction;
    const auto gsfBackward = gsfForward.invert();

    // Check if the start parameters are on the start surface
    auto intersectionStatusStartSurface =
        sParameters.referenceSurface()
            .intersect(GeometryContext{},
                       sParameters.position(GeometryContext{}),
                       sParameters.unitDirection(), true)
            .intersection.status;

    if (intersectionStatusStartSurface != Intersection3D::Status::onSurface) {
      ACTS_ERROR(
          "Surface intersection of start parameters with bound-check failed");
      return GsfError::StartParametersNotOnStartSurface;
    }

    // To be able to find measurements later, we put them into a map
    // We need to copy input SourceLinks anyways, so the map can own them.
    ACTS_VERBOSE("Preparing " << std::distance(begin, end)
                              << " input measurements");
    std::map<GeometryIdentifier, SourceLink> inputMeasurements;
    for (auto it = begin; it != end; ++it) {
      SourceLink sl = *it;
      auto geoId = sl.geometryId();
      inputMeasurements.emplace(geoId, std::move(sl));
    }

    ACTS_VERBOSE(
        "Gsf: Final measuerement map size: " << inputMeasurements.size());

    if (sParameters.covariance() == std::nullopt) {
      return GsfError::StartParametersHaveNoCovariance;
    }

    /////////////////
    // Forward pass
    /////////////////
    ACTS_VERBOSE("+-----------------------------+");
    ACTS_VERBOSE("| Gsf: Do forward propagation |");
    ACTS_VERBOSE("+-----------------------------+");

    auto fwdResult = [&]() {
      auto fwdPropOptions = fwdPropInitializer(options);

      // Catch the actor and set the measurements
      auto& actor = fwdPropOptions.actionList.template get<GsfActor>();
      actor.setOptions(options);
      actor.m_cfg.inputMeasurements = &inputMeasurements;
      actor.m_cfg.numberMeasurements = inputMeasurements.size();
      actor.m_cfg.inReversePass = false;
      actor.m_cfg.logger = m_actorLogger.get();

      fwdPropOptions.direction = gsfForward;

      // If necessary convert to MultiComponentBoundTrackParameters
      using IsMultiParameters =
          detail::IsMultiComponentBoundParameters<start_parameters_t>;

      typename propagator_t::template action_list_t_result_t<
          CurvilinearTrackParameters, decltype(fwdPropOptions.actionList)>
          inputResult;

      auto& r = inputResult.template get<detail::GsfResult<traj_t>>();

      r.fittedStates = &trackContainer.trackStateContainer();

      // This allows the initialization with single- and multicomponent start
      // parameters
      if constexpr (not IsMultiParameters::value) {
        using Charge = typename IsMultiParameters::Charge;

        MultiComponentBoundTrackParameters<Charge> params(
            sParameters.referenceSurface().getSharedPtr(),
            sParameters.parameters(), sParameters.covariance());

        return m_propagator.propagate(params, fwdPropOptions,
                                      std::move(inputResult));
      } else {
        return m_propagator.propagate(sParameters, fwdPropOptions,
                                      std::move(inputResult));
      }
    }();

    if (!fwdResult.ok()) {
      return return_error_or_abort(fwdResult.error());
    }

    auto& fwdGsfResult = fwdResult->template get<detail::GsfResult<traj_t>>();

    if (!fwdGsfResult.result.ok()) {
      return return_error_or_abort(fwdGsfResult.result.error());
    }

    if (fwdGsfResult.measurementStates == 0) {
      return return_error_or_abort(GsfError::NoMeasurementStatesCreatedForward);
    }

    ACTS_VERBOSE("Finished forward propagation");
    ACTS_VERBOSE("- visited surfaces: " << fwdGsfResult.visitedSurfaces.size());
    ACTS_VERBOSE("- processed states: " << fwdGsfResult.processedStates);
    ACTS_VERBOSE("- measuerement states: " << fwdGsfResult.measurementStates);

    //////////////////
    // Backward pass
    //////////////////
    ACTS_VERBOSE("+------------------------------+");
    ACTS_VERBOSE("| Gsf: Do backward propagation |");
    ACTS_VERBOSE("+------------------------------+");

    auto bwdResult = [&]() {
      auto bwdPropOptions = bwdPropInitializer(options);

      auto& actor = bwdPropOptions.actionList.template get<GsfActor>();
      actor.setOptions(options);
      actor.m_cfg.inputMeasurements = &inputMeasurements;
      actor.m_cfg.inReversePass = true;
      actor.m_cfg.logger = m_actorLogger.get();
      actor.setOptions(options);

      bwdPropOptions.direction = gsfBackward;

      const Surface& target = options.referenceSurface
                                  ? *options.referenceSurface
                                  : sParameters.referenceSurface();

      using PM = TrackStatePropMask;

      typename propagator_t::template action_list_t_result_t<
          BoundTrackParameters, decltype(bwdPropOptions.actionList)>
          inputResult;

      // Unfortunately we must construct the result type here to be able to
      // return an error code
      using ResultType =
          decltype(m_propagator.template propagate<
                   MultiComponentBoundTrackParameters<SinglyCharged>,
                   decltype(bwdPropOptions), MultiStepperSurfaceReached>(
              std::declval<MultiComponentBoundTrackParameters<SinglyCharged>>(),
              std::declval<Acts::Surface&>(),
              std::declval<decltype(bwdPropOptions)>(),
              std::declval<decltype(inputResult)>()));

      auto& r = inputResult.template get<detail::GsfResult<traj_t>>();

      r.fittedStates = &trackContainer.trackStateContainer();

      assert(
          (fwdGsfResult.lastMeasurementTip != MultiTrajectoryTraits::kInvalid &&
           "tip is invalid"));

      auto proxy =
          r.fittedStates->getTrackState(fwdGsfResult.lastMeasurementTip);
      proxy.shareFrom(TrackStatePropMask::Filtered,
                      TrackStatePropMask::Smoothed);

      r.currentTip = fwdGsfResult.lastMeasurementTip;
      r.visitedSurfaces.push_back(&proxy.referenceSurface());
      r.surfacesVisitedBwdAgain.push_back(&proxy.referenceSurface());
      r.measurementStates++;
      r.processedStates++;

      const auto& params = *fwdGsfResult.lastMeasurementState;

      return m_propagator.template propagate<std::decay_t<decltype(params)>,
                                             decltype(bwdPropOptions),
                                             MultiStepperSurfaceReached>(
          params, target, bwdPropOptions, std::move(inputResult));
    }();

    if (!bwdResult.ok()) {
      return return_error_or_abort(bwdResult.error());
    }

    auto& bwdGsfResult = bwdResult->template get<detail::GsfResult<traj_t>>();

    if (!bwdGsfResult.result.ok()) {
      return return_error_or_abort(bwdGsfResult.result.error());
    }

    if (bwdGsfResult.measurementStates == 0) {
      return return_error_or_abort(
          GsfError::NoMeasurementStatesCreatedBackward);
    }

    ////////////////////////////////////
    // Create Kalman Result
    ////////////////////////////////////
    ACTS_VERBOSE("Gsf - States summary:");
    ACTS_VERBOSE("- Fwd measurement states: " << fwdGsfResult.measurementStates
                                              << ", holes: "
                                              << fwdGsfResult.measurementHoles);
    ACTS_VERBOSE("- Bwd measurement states: " << bwdGsfResult.measurementStates
                                              << ", holes: "
                                              << bwdGsfResult.measurementHoles);

    // TODO should this be warning level? it happens quite often... Investigate!
    if (bwdGsfResult.measurementStates != fwdGsfResult.measurementStates) {
      ACTS_DEBUG("Fwd and bwd measuerement states do not match");
    }

    // Go through the states and assign outliers / unset smoothed if surface not
    // passed in backward pass
    const auto& foundBwd = bwdGsfResult.surfacesVisitedBwdAgain;
    std::size_t measurementStatesFinal = 0;

    for (auto state :
         fwdGsfResult.fittedStates->trackStateRange(fwdGsfResult.currentTip)) {
      const bool found = std::find(foundBwd.begin(), foundBwd.end(),
                                   &state.referenceSurface()) != foundBwd.end();
      if (not found && state.typeFlags().test(MeasurementFlag)) {
        state.typeFlags().set(OutlierFlag);
        state.typeFlags().reset(MeasurementFlag);
        state.unset(TrackStatePropMask::Smoothed);
      }

      measurementStatesFinal +=
          static_cast<std::size_t>(state.typeFlags().test(MeasurementFlag));
    }

    if (measurementStatesFinal == 0) {
      return return_error_or_abort(GsfError::NoMeasurementStatesCreatedFinal);
    }

    auto track = trackContainer.getTrack(trackContainer.addTrack());
    track.tipIndex() = fwdGsfResult.lastMeasurementTip;

    if (options.referenceSurface) {
      const auto& params = *bwdResult->endParameters;
      track.parameters() = params.parameters();
      track.covariance() = params.covariance().value();
      track.setReferenceSurface(params.referenceSurface().getSharedPtr());

      if (trackContainer.hasColumn(
              hashString(GsfConstants::kFinalMultiComponentStateColumn))) {
        ACTS_DEBUG("Add final multi-component state to track")
        const auto& fsr = bwdResult->template get<
            Acts::detail::FinalStateCollector::result_type>();
        track.template component<GsfConstants::FinalMultiComponentState>(
            GsfConstants::kFinalMultiComponentStateColumn) = fsr.pars;
      }
    }

    calculateTrackQuantities(track);

    track.nMeasurements() = measurementStatesFinal;
    track.nHoles() = fwdGsfResult.measurementHoles;

    return track;
  }
};

}  // namespace Experimental
}  // namespace Acts

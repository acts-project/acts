// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/MultiStepperAborters.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/BetheHeitlerApprox.hpp"
#include "Acts/TrackFitting/detail/GsfActor.hpp"

#include <fstream>

namespace Acts {

namespace detail {

/// Type trait to identify if a type is a MultiComponentBoundTrackParameters and
/// to inspect its charge representation if not TODO this probably gives an ugly
/// error message if detectCharge does not compile
template <typename T>
struct IsMultiComponentBoundParameters : public std::false_type {
  template <template <class> class U, class V>
  static auto detectCharge(const U<V>&) {
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
///
/// @note This GSF implementation tries to be as compatible to the KalmanFitter
/// as possible. However, there are certain differences at the moment:
/// * There is always a backward pass during fitting.
/// * There are only measurement states in the result
/// * Passed-again-surfaces is always empty at the moment
/// * Probably some more differences which I don't think of at the moment.
template <typename propagator_t, typename traj_t,
          typename bethe_heitler_approx_t = detail::BetheHeitlerApprox<6, 5>>
struct GaussianSumFitter {
  GaussianSumFitter(propagator_t&& propagator,
                    bethe_heitler_approx_t&& bha = bethe_heitler_approx_t(
                        detail::bh_cdf_cmps6_order5_data))
      : m_propagator(std::move(propagator)), m_bethe_heitler_approx(bha) {}

  /// The propagator instance used by the fit function
  propagator_t m_propagator;

  /// The fitter holds the instance of the bethe heitler approx
  bethe_heitler_approx_t m_bethe_heitler_approx;

  /// The navigator type
  using GsfNavigator = typename propagator_t::Navigator;

  /// The actor type
  using GsfActor = detail::GsfActor<bethe_heitler_approx_t, traj_t>;

  /// @brief The fit function for the Direct navigator
  template <typename source_link_it_t, typename start_parameters_t>
  Acts::Result<Acts::KalmanFitterResult<traj_t>> fit(
      source_link_it_t begin, source_link_it_t end,
      const start_parameters_t& sParameters, const GsfOptions<traj_t>& options,
      const std::vector<const Surface*>& sSequence,
      std::shared_ptr<traj_t> trajectory = {}) const {
    // Check if we have the correct navigator
    static_assert(
        std::is_same_v<DirectNavigator, typename propagator_t::Navigator>);

    // Initialize the forward propagation with the DirectNavigator
    auto fwdPropInitializer = [&sSequence, this](const auto& opts,
                                                 const auto& logger) {
      using Actors = ActionList<GsfActor, DirectNavigator::Initializer>;
      using Aborters = AbortList<>;

      PropagatorOptions<Actors, Aborters> propOptions(
          opts.geoContext, opts.magFieldContext, logger);

      propOptions.setPlainOptions(opts.propagatorPlainOptions);

      propOptions.actionList.template get<DirectNavigator::Initializer>()
          .navSurfaces = sSequence;
      propOptions.actionList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_bethe_heitler_approx;

      return propOptions;
    };

    // Initialize the backward propagation with the DirectNavigator
    auto bwdPropInitializer = [&sSequence, this](const auto& opts,
                                                 const auto& logger) {
      using Actors = ActionList<GsfActor, DirectNavigator::Initializer>;
      using Aborters = AbortList<>;

      std::vector<const Surface*> backwardSequence(
          std::next(sSequence.rbegin()), sSequence.rend());
      backwardSequence.push_back(opts.referenceSurface);

      PropagatorOptions<Actors, Aborters> propOptions(
          opts.geoContext, opts.magFieldContext, logger);

      propOptions.setPlainOptions(opts.propagatorPlainOptions);

      propOptions.actionList.template get<DirectNavigator::Initializer>()
          .navSurfaces = std::move(backwardSequence);
      propOptions.actionList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_bethe_heitler_approx;

      return propOptions;
    };

    if (!trajectory) {
      trajectory = std::make_shared<traj_t>();
    }

    return fit_impl(begin, end, sParameters, options, fwdPropInitializer,
                    bwdPropInitializer, std::move(trajectory));
  }

  /// @brief The fit function for the standard navigator
  template <typename source_link_it_t, typename start_parameters_t>
  Acts::Result<Acts::KalmanFitterResult<traj_t>> fit(
      source_link_it_t begin, source_link_it_t end,
      const start_parameters_t& sParameters, const GsfOptions<traj_t>& options,
      std::shared_ptr<traj_t> trajectory = {}) const {
    // Check if we have the correct navigator
    static_assert(std::is_same_v<Navigator, typename propagator_t::Navigator>);

    // Initialize the forward propagation with the DirectNavigator
    auto fwdPropInitializer = [this](const auto& opts, const auto& logger) {
      using Actors = ActionList<GsfActor>;
      using Aborters = AbortList<EndOfWorldReached>;

      PropagatorOptions<Actors, Aborters> propOptions(
          opts.geoContext, opts.magFieldContext, logger);
      propOptions.setPlainOptions(opts.propagatorPlainOptions);
      propOptions.actionList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_bethe_heitler_approx;

      return propOptions;
    };

    // Initialize the backward propagation with the DirectNavigator
    auto bwdPropInitializer = [this](const auto& opts, const auto& logger) {
      using Actors = ActionList<GsfActor>;
      using Aborters = AbortList<EndOfWorldReached>;

      PropagatorOptions<Actors, Aborters> propOptions(
          opts.geoContext, opts.magFieldContext, logger);

      propOptions.setPlainOptions(opts.propagatorPlainOptions);

      propOptions.actionList.template get<GsfActor>()
          .m_cfg.bethe_heitler_approx = &m_bethe_heitler_approx;

      return propOptions;
    };

    if (!trajectory) {
      trajectory = std::make_shared<traj_t>();
    }
    return fit_impl(begin, end, sParameters, options, fwdPropInitializer,
                    bwdPropInitializer, std::move(trajectory));
  }

  /// The generic implementation of the fit function.
  /// TODO check what this function does with the referenceSurface is e.g. the
  /// first measuerementSurface
  template <typename source_link_it_t, typename start_parameters_t,
            typename fwd_prop_initializer_t, typename bwd_prop_initializer_t>
  Acts::Result<Acts::KalmanFitterResult<traj_t>> fit_impl(
      source_link_it_t begin, source_link_it_t end,
      const start_parameters_t& sParameters, const GsfOptions<traj_t>& options,
      const fwd_prop_initializer_t& fwdPropInitializer,
      const bwd_prop_initializer_t& bwdPropInitializer,
      std::shared_ptr<traj_t> trajectory) const {
    // return or abort utility
    auto return_error_or_abort = [&](auto error) {
      if (options.abortOnError) {
        std::abort();
      }
      return error;
    };

    // The logger
    const auto& logger = options.logger;

    // Define directions based on input propagation direction. This way we can
    // refer to 'forward' and 'backward' regardless of the actual direction.
    const auto gsfForward = options.propagatorPlainOptions.direction;
    const auto gsfBackward = static_cast<NavigationDirection>(-1 * gsfForward);

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
    std::map<GeometryIdentifier, std::reference_wrapper<const SourceLink>>
        inputMeasurements;
    for (auto it = begin; it != end; ++it) {
      const SourceLink& sl = *it;
      inputMeasurements.emplace(sl.geometryId(), *it);
    }

    ACTS_VERBOSE(
        "Gsf: Final measuerement map size: " << inputMeasurements.size());
    throw_assert(sParameters.covariance() != std::nullopt,
                 "we need a covariance here...");

    /////////////////
    // Forward pass
    /////////////////
    ACTS_VERBOSE("+-----------------------------+");
    ACTS_VERBOSE("| Gsf: Do forward propagation |");
    ACTS_VERBOSE("+-----------------------------+");

    auto fwdResult = [&]() {
      auto fwdPropOptions = fwdPropInitializer(options, logger);

      // Catch the actor and set the measurements
      auto& actor = fwdPropOptions.actionList.template get<GsfActor>();
      actor.m_cfg.inputMeasurements = inputMeasurements;
      actor.m_cfg.maxComponents = options.maxComponents;
      actor.m_cfg.extensions = options.extensions;
      actor.m_cfg.abortOnError = options.abortOnError;
      actor.m_cfg.disableAllMaterialHandling =
          options.disableAllMaterialHandling;
      actor.m_cfg.numberMeasurements = inputMeasurements.size();

      fwdPropOptions.direction = gsfForward;

      // If necessary convert to MultiComponentBoundTrackParameters
      using IsMultiParameters =
          detail::IsMultiComponentBoundParameters<start_parameters_t>;

      typename propagator_t::template action_list_t_result_t<
          CurvilinearTrackParameters, decltype(fwdPropOptions.actionList)>
          inputResult;

      auto& r = inputResult.template get<detail::GsfResult<traj_t>>();

      r.fittedStates = trajectory;

      // This allows the initialization with single- and multicomponent start
      // parameters
      if constexpr (not IsMultiParameters::value) {
        using Charge = typename IsMultiParameters::Charge;

        r.parentTips.resize(1, MultiTrajectoryTraits::kInvalid);

        MultiComponentBoundTrackParameters<Charge> params(
            sParameters.referenceSurface().getSharedPtr(),
            sParameters.parameters(), sParameters.covariance());

        return m_propagator.propagate(params, fwdPropOptions,
                                      std::move(inputResult));
      } else {
        r.parentTips.resize(sParameters.components().size(),
                            MultiTrajectoryTraits::kInvalid);

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

    if (fwdGsfResult.processedStates == 0) {
      return return_error_or_abort(GsfError::NoStatesCreated);
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
      auto bwdPropOptions = bwdPropInitializer(options, logger);

      auto& actor = bwdPropOptions.actionList.template get<GsfActor>();
      actor.m_cfg.inputMeasurements = inputMeasurements;
      actor.m_cfg.maxComponents = options.maxComponents;
      actor.m_cfg.abortOnError = options.abortOnError;
      actor.m_cfg.disableAllMaterialHandling =
          options.disableAllMaterialHandling;
      actor.m_cfg.extensions = options.extensions;

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

      r.fittedStates = trajectory;

      // We take the last measurement state (filtered) from the forward result
      // as the first measurement state in the backward result (predicted and
      // filtered), so we can skip the Kalman update on the first surface as
      // this would be redundant. We combine this with the construction of the
      // propagation start parameters to ensure they are consistent.
      std::vector<std::tuple<double, BoundVector, BoundSymMatrix>> cmps;
      std::shared_ptr<const Surface> surface;

      for (const auto idx : fwdGsfResult.lastMeasurementTips) {
        // TODO This should not happen, but very rarely does. Maybe investigate
        // later
        if (fwdGsfResult.weightsOfStates.at(idx) == 0) {
          continue;
        }

        r.currentTips.push_back(
            r.fittedStates->addTrackState(TrackStatePropMask::All));

        auto proxy = r.fittedStates->getTrackState(r.currentTips.back());
        proxy.copyFrom(fwdGsfResult.fittedStates->getTrackState(idx));
        r.weightsOfStates[r.currentTips.back()] =
            fwdGsfResult.weightsOfStates.at(idx);

        proxy.shareFrom(proxy, PM::Filtered, PM::Predicted);

        // Avoid accessing the surface for every component, since it should be
        // the same
        if (not surface) {
          surface = proxy.referenceSurface().getSharedPtr();
        }

        cmps.push_back({fwdGsfResult.weightsOfStates.at(idx), proxy.filtered(),
                        proxy.filteredCovariance()});
      }

      if (cmps.empty()) {
        return ResultType{GsfError::NoComponentCreated};
      }

      r.visitedSurfaces.insert(surface->geometryId());
      r.parentTips = r.currentTips;
      r.measurementStates++;
      r.processedStates++;

      const auto params =
          MultiComponentBoundTrackParameters<SinglyCharged>(surface, cmps);

      return m_propagator
          .template propagate<decltype(params), decltype(bwdPropOptions),
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

    if (bwdGsfResult.processedStates == 0) {
      return return_error_or_abort(GsfError::NoStatesCreated);
    }

    ////////////////////////////////////
    // Smooth and create Kalman Result
    ////////////////////////////////////
    ACTS_VERBOSE("Gsf: Do smoothing");
    ACTS_VERBOSE("- Fwd measurement states: " << fwdGsfResult.measurementStates
                                              << ", holes: "
                                              << fwdGsfResult.measurementHoles);
    ACTS_VERBOSE("- Bwd measurement states: " << bwdGsfResult.measurementStates
                                              << ", holes: "
                                              << bwdGsfResult.measurementHoles);

    auto smoothResult = detail::smoothAndCombineTrajectories<traj_t, true>(
        *fwdGsfResult.fittedStates, fwdGsfResult.currentTips,
        fwdGsfResult.weightsOfStates, *bwdGsfResult.fittedStates,
        bwdGsfResult.currentTips, bwdGsfResult.weightsOfStates, logger);

    // Cannot use structured binding since they cannot be captured in lambda
    auto& kalmanResult = std::get<0>(smoothResult);

    // Some test
    if (std::get<1>(smoothResult).empty()) {
      return return_error_or_abort(GsfError::NoStatesCreated);
    }

    // Compute the missed active surfaces as the union of the forward and
    // backward pass missed active surfaces
    // TODO this is quite expencive computationally, maybe just use from fwd?
    {
      auto fwdActSurf = fwdGsfResult.missedActiveSurfaces;
      std::sort(fwdActSurf.begin(), fwdActSurf.end());

      auto bwdActSurf = bwdGsfResult.missedActiveSurfaces;
      std::sort(bwdActSurf.begin(), bwdActSurf.end());

      std::vector<const Surface*> missedActiveSurfaces;
      std::set_union(fwdActSurf.begin(), fwdActSurf.end(), bwdActSurf.begin(),
                     bwdActSurf.end(),
                     std::back_inserter(missedActiveSurfaces));

      kalmanResult.missedActiveSurfaces = missedActiveSurfaces;
    }

    //////////////////////////////////////////////////////////////////
    // Propagate back to reference surface with smoothed parameters //
    //////////////////////////////////////////////////////////////////
    if (options.referenceSurface) {
      ACTS_VERBOSE("+-----------------------------------------------+");
      ACTS_VERBOSE("| Gsf: Do propagation back to reference surface |");
      ACTS_VERBOSE("+-----------------------------------------------+");
      auto lastResult = [&]() -> Result<std::unique_ptr<BoundTrackParameters>> {
        const auto& [surface, lastSmoothedState] =
            std::get<1>(smoothResult).front();

        throw_assert(
            detail::weightsAreNormalized(
                lastSmoothedState,
                [](const auto& tuple) { return std::get<double>(tuple); }),
            "");

        const MultiComponentBoundTrackParameters<SinglyCharged> params(
            surface->getSharedPtr(), lastSmoothedState);

        auto lastPropOptions = bwdPropInitializer(options, logger);

        auto& actor = lastPropOptions.actionList.template get<GsfActor>();
        actor.m_cfg.maxComponents = options.maxComponents;
        actor.m_cfg.abortOnError = options.abortOnError;
        actor.m_cfg.disableAllMaterialHandling =
            options.disableAllMaterialHandling;

        lastPropOptions.direction = gsfBackward;

        typename propagator_t::template action_list_t_result_t<
            BoundTrackParameters, decltype(lastPropOptions.actionList)>
            inputResult;

        auto& r = inputResult.template get<detail::GsfResult<traj_t>>();

        r.fittedStates = trajectory;
        r.parentTips.resize(params.components().size(),
                            MultiTrajectoryTraits::kInvalid);

        // Add the initial surface to the list of already visited surfaces, so
        // that the material effects are not applied twice
        r.visitedSurfaces.insert(surface->geometryId());

        auto result =
            m_propagator
                .template propagate<decltype(params), decltype(lastPropOptions),
                                    MultiStepperSurfaceReached>(
                    params, *options.referenceSurface, lastPropOptions,
                    std::move(inputResult));

        if (!result.ok()) {
          return result.error();
        } else {
          return std::move(result->endParameters);
        }
      }();

      if (!lastResult.ok()) {
        return return_error_or_abort(lastResult.error());
      }

      kalmanResult.fittedParameters = **lastResult;
    }

    return kalmanResult;
  }
};

}  // namespace Experimental
}  // namespace Acts

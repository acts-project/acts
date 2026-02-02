// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"
#include "Acts/TrackFitting/detail/GsfComponentMerging.hpp"
#include "Acts/TrackFitting/detail/GsfUtils.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <map>

namespace Acts::detail {

template <typename traj_t>
struct GsfResult {
  /// The multi-trajectory which stores the graph of components
  traj_t* fittedStates{nullptr};

  /// The current top index of the MultiTrajectory
  TrackIndexType currentTip = kTrackIndexInvalid;

  /// The last tip referring to a measurement state in the MultiTrajectory
  TrackIndexType lastMeasurementTip = kTrackIndexInvalid;

  /// The last multi-component measurement state. Used to initialize the
  /// backward pass.
  std::vector<std::tuple<double, BoundVector, BoundMatrix>>
      lastMeasurementComponents;

  /// The last measurement surface. Used to initialize the backward pass.
  const Acts::Surface* lastMeasurementSurface = nullptr;

  /// Some counting
  std::size_t measurementStates = 0;
  std::size_t measurementHoles = 0;
  std::size_t processedStates = 0;

  std::vector<const Surface*> visitedSurfaces;
  std::vector<const Surface*> surfacesVisitedBwdAgain;

  /// Statistics about material encounterings
  Updatable<std::size_t> nInvalidBetheHeitler;
  Updatable<double> maxPathXOverX0;
  Updatable<double> sumPathXOverX0;

  // Internal: bethe heitler approximation component cache
  std::vector<BetheHeitlerApprox::Component> betheHeitlerCache;

  // Internal: component cache to avoid reallocation
  std::vector<GsfComponent> componentCache;
};

/// The actor carrying out the GSF algorithm
template <typename traj_t>
struct GsfActor {
  /// Enforce default construction
  GsfActor() = default;

  using ComponentCache = GsfComponent;

  /// Broadcast the result_type
  using result_type = GsfResult<traj_t>;

  // Actor configuration
  struct Config {
    /// Maximum number of components which the GSF should handle
    std::size_t maxComponents = 16;

    /// Input measurements
    const std::map<GeometryIdentifier, SourceLink>* inputMeasurements = nullptr;

    /// Bethe Heitler Approximator pointer. The fitter holds the approximator
    /// instance TODO if we somehow could initialize a reference here...
    const BetheHeitlerApprox* bethe_heitler_approx = nullptr;

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// When to discard components
    double weightCutoff = 1.0e-4;

    /// When this option is enabled, material information on all surfaces is
    /// ignored. This disables the component convolution as well as the handling
    /// of energy. This may be useful for debugging.
    bool disableAllMaterialHandling = false;

    /// Whether to abort immediately when an error occurs
    bool abortOnError = false;

    /// We can stop the propagation if we reach this number of measurement
    /// states
    std::optional<std::size_t> numberMeasurements;

    /// The extensions
    GsfExtensions<traj_t> extensions;

    /// Whether we are in the reverse pass or not. This is more reliable than
    /// checking the navigation direction, because in principle the fitter can
    /// be started backwards in the first pass
    bool inReversePass = false;

    /// How to reduce the states that are stored in the multi trajectory
    ComponentMergeMethod mergeMethod = ComponentMergeMethod::eMaxWeight;

    const Logger* logger{nullptr};

    /// Calibration context for the fit
    const CalibrationContext* calibrationContext{nullptr};

  } m_cfg;

  const Logger& logger() const { return *m_cfg.logger; }

  struct TemporaryStates {
    traj_t traj;
    std::vector<TrackIndexType> tips;
    std::map<TrackIndexType, double> weights;
  };

  using FiltProjector = MultiTrajectoryProjector<StatesType::eFiltered, traj_t>;

  /// @brief GSF actor operation
  ///
  /// @tparam propagator_state_t is the type of Propagator state
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param state is the mutable propagator state object
  /// @param stepper The stepper in use
  /// @param result is the mutable result state object
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  Result<void> act(propagator_state_t& state, const stepper_t& stepper,
                   const navigator_t& navigator, result_type& result,
                   const Logger& /*logger*/) const {
    assert(result.fittedStates && "No MultiTrajectory set");

    // Prints some VERBOSE things and performs some asserts. Can be removed
    // without change of behaviour
    const detail::ScopedGsfInfoPrinterAndChecker printer(state, stepper,
                                                         navigator, logger());

    // We only need to do something if we are on a surface
    if (!navigator.currentSurface(state.navigation)) {
      return Result<void>::success();
    }

    const auto& surface = *navigator.currentSurface(state.navigation);
    ACTS_VERBOSE("Step is at surface " << surface.geometryId());

    // All components must be normalized at the beginning here, otherwise the
    // stepper misbehaves
    [[maybe_unused]] auto stepperComponents =
        stepper.constComponentIterable(state.stepping);
    assert(detail::weightsAreNormalized(
        stepperComponents, [](const auto& cmp) { return cmp.weight(); }));

    // All components must have status "on surface". It is however possible,
    // that currentSurface is nullptr and all components are "on surface" (e.g.,
    // for surfaces excluded from the navigation)
    using Status [[maybe_unused]] = IntersectionStatus;
    assert(std::all_of(
        stepperComponents.begin(), stepperComponents.end(),
        [](const auto& cmp) { return cmp.status() == Status::onSurface; }));

    // Early return if we already were on this surface TODO why is this
    // necessary
    const bool visited = rangeContainsValue(result.visitedSurfaces, &surface);

    if (visited) {
      ACTS_VERBOSE("Already visited surface, return");
      return Result<void>::success();
    }

    result.visitedSurfaces.push_back(&surface);

    // Check what we have on this surface
    const auto foundSourceLink =
        m_cfg.inputMeasurements->find(surface.geometryId());
    const bool haveMaterial =
        navigator.currentSurface(state.navigation)->surfaceMaterial() &&
        !m_cfg.disableAllMaterialHandling;
    const bool haveMeasurement =
        foundSourceLink != m_cfg.inputMeasurements->end();

    ACTS_VERBOSE(std::boolalpha << "haveMaterial " << haveMaterial
                                << ", haveMeasurement: " << haveMeasurement);

    ////////////////////////
    // The Core Algorithm
    ////////////////////////

    // Early return if nothing happens
    if (!haveMaterial && !haveMeasurement) {
      // No hole before first measurement
      if (result.processedStates > 0 && surface.isSensitive()) {
        TemporaryStates tmpStates;
        noMeasurementUpdate(state, stepper, navigator, result, tmpStates, true);
      }
      return Result<void>::success();
    }

    // Update the counters. Note that this should be done before potential
    // material interactions, because if this is our last measurement this would
    // not influence the fit anymore.
    if (haveMeasurement) {
      result.maxPathXOverX0.update();
      result.sumPathXOverX0.update();
      result.nInvalidBetheHeitler.update();
    }

    for (auto cmp : stepper.componentIterable(state.stepping)) {
      cmp.singleStepper(stepper).transportCovarianceToBound(cmp.state(),
                                                            surface);
    }

    if (haveMaterial) {
      if (haveMeasurement) {
        applyMultipleScattering(
            state, stepper, navigator,
            detail::determineMaterialUpdateMode(state, navigator,
                                                MaterialUpdateMode::PreUpdate));
      } else {
        applyMultipleScattering(
            state, stepper, navigator,
            detail::determineMaterialUpdateMode(
                state, navigator, MaterialUpdateMode::FullUpdate));
      }
    }

    // We do not need the component cache here, we can just update our stepper
    // state with the filtered components.
    // NOTE because of early return before we know that we have a measurement
    if (!haveMaterial) {
      TemporaryStates tmpStates;

      auto res = kalmanUpdate(state, stepper, navigator, result, tmpStates,
                              foundSourceLink->second);

      if (!res.ok()) {
        if (m_cfg.abortOnError) {
          std::abort();
        }
        return res.error();
      }

      updateStepper(state, stepper, tmpStates);
    }
    // We have material, we thus need a component cache since we will
    // convolute the components and later reduce them again before updating
    // the stepper
    else {
      TemporaryStates tmpStates;
      Result<void> res;

      if (haveMeasurement) {
        res = kalmanUpdate(state, stepper, navigator, result, tmpStates,
                           foundSourceLink->second);
      } else {
        res = noMeasurementUpdate(state, stepper, navigator, result, tmpStates,
                                  false);
      }

      if (!res.ok()) {
        if (m_cfg.abortOnError) {
          std::abort();
        }
        return res.error();
      }

      // Reuse memory over all calls to the Actor in a single propagation
      std::vector<ComponentCache>& componentCache = result.componentCache;
      componentCache.clear();

      convoluteComponents(state, stepper, navigator, tmpStates, componentCache,
                          result);

      if (componentCache.empty()) {
        ACTS_WARNING(
            "No components left after applying energy loss. "
            "Is the weight cutoff "
            << m_cfg.weightCutoff << " too high?");
        ACTS_WARNING("Return to propagator without applying energy loss");
        return Result<void>::success();
      }

      // reduce component number
      const auto finalCmpNumber = std::min(
          static_cast<std::size_t>(stepper.maxComponents), m_cfg.maxComponents);
      m_cfg.extensions.mixtureReducer(componentCache, finalCmpNumber, surface);

      removeLowWeightComponents(componentCache);

      updateStepper(state, stepper, navigator, componentCache);
    }

    // If we have only done preUpdate before, now do postUpdate
    if (haveMaterial && haveMeasurement) {
      applyMultipleScattering(
          state, stepper, navigator,
          detail::determineMaterialUpdateMode(state, navigator,
                                              MaterialUpdateMode::PostUpdate));
    }

    return Result<void>::success();
  }

  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool checkAbort(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                  const navigator_t& /*navigator*/, const result_type& result,
                  const Logger& /*logger*/) const {
    if (m_cfg.numberMeasurements &&
        result.measurementStates == m_cfg.numberMeasurements) {
      ACTS_VERBOSE("Stop navigation because all measurements are found");
      return true;
    }

    return false;
  }

  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  void convoluteComponents(propagator_state_t& state, const stepper_t& stepper,
                           const navigator_t& navigator,
                           const TemporaryStates& tmpStates,
                           std::vector<ComponentCache>& componentCache,
                           result_type& result) const {
    double pathXOverX0 = 0.0;
    auto cmps = stepper.componentIterable(state.stepping);
    for (auto [idx, cmp] : zip(tmpStates.tips, cmps)) {
      auto proxy = tmpStates.traj.getTrackState(idx);

      BoundTrackParameters bound(proxy.referenceSurface().getSharedPtr(),
                                 proxy.filtered(), proxy.filteredCovariance(),
                                 stepper.particleHypothesis(state.stepping));

      pathXOverX0 +=
          applyBetheHeitler(state, navigator, bound, tmpStates.weights.at(idx),
                            componentCache, result);
    }

    // Store average material seen by the components
    // Should not be too broadly distributed
    result.sumPathXOverX0.tmp() += pathXOverX0 / tmpStates.tips.size();
  }

  template <typename propagator_state_t, typename navigator_t>
  double applyBetheHeitler(const propagator_state_t& state,
                           const navigator_t& navigator,
                           const BoundTrackParameters& old_bound,
                           const double old_weight,
                           std::vector<ComponentCache>& componentCache,
                           result_type& result) const {
    const auto& surface = *navigator.currentSurface(state.navigation);
    const auto p_prev = old_bound.absoluteMomentum();
    const auto& particleHypothesis = old_bound.particleHypothesis();

    // Evaluate material slab
    auto slab = surface.surfaceMaterial()->materialSlab(
        old_bound.position(state.geoContext), state.options.direction,
        MaterialUpdateMode::FullUpdate);

    const auto pathCorrection = surface.pathCorrection(
        state.geoContext, old_bound.position(state.geoContext),
        old_bound.direction());
    slab.scaleThickness(pathCorrection);

    const double pathXOverX0 = slab.thicknessInX0();
    result.maxPathXOverX0.tmp() =
        std::max(result.maxPathXOverX0.tmp(), pathXOverX0);

    // Emit a warning if the approximation is not valid for this x/x0
    if (!m_cfg.bethe_heitler_approx->validXOverX0(pathXOverX0)) {
      ++result.nInvalidBetheHeitler.tmp();
      ACTS_DEBUG(
          "Bethe-Heitler approximation encountered invalid value for x/x0="
          << pathXOverX0 << " at surface " << surface.geometryId());
    }

    // Get the mixture
    result.betheHeitlerCache.resize(
        m_cfg.bethe_heitler_approx->maxComponents());
    const auto mixture = m_cfg.bethe_heitler_approx->mixture(
        pathXOverX0, result.betheHeitlerCache);

    // Create all possible new components
    for (const auto& gaussian : mixture) {
      // Here we combine the new child weight with the parent weight.
      // However, this must be later re-adjusted
      const auto new_weight = gaussian.weight * old_weight;

      if (new_weight < m_cfg.weightCutoff) {
        ACTS_VERBOSE("Skip component with weight " << new_weight);
        continue;
      }

      if (gaussian.mean < 1.e-8) {
        ACTS_WARNING("Skip component with gaussian " << gaussian.mean << " +- "
                                                     << gaussian.var);
        continue;
      }

      // compute delta p from mixture and update parameters
      auto new_pars = old_bound.parameters();

      const auto delta_p = [&]() {
        if (state.options.direction == Direction::Forward()) {
          return p_prev * (gaussian.mean - 1.);
        } else {
          return p_prev * (1. / gaussian.mean - 1.);
        }
      }();

      assert(p_prev + delta_p > 0. && "new momentum must be > 0");
      new_pars[eBoundQOverP] =
          particleHypothesis.qOverP(p_prev + delta_p, old_bound.charge());

      // compute inverse variance of p from mixture and update covariance
      auto new_cov = old_bound.covariance().value();

      const auto varInvP = [&]() {
        if (state.options.direction == Direction::Forward()) {
          const auto f = 1. / (p_prev * gaussian.mean);
          return f * f * gaussian.var;
        } else {
          return gaussian.var / (p_prev * p_prev);
        }
      }();

      new_cov(eBoundQOverP, eBoundQOverP) += varInvP;
      assert(std::isfinite(new_cov(eBoundQOverP, eBoundQOverP)) &&
             "new cov not finite");

      // Set the remaining things and push to vector
      componentCache.push_back({new_weight, new_pars, new_cov});
    }

    return pathXOverX0;
  }

  /// Remove components with low weights and renormalize from the component
  /// cache
  /// TODO This function does not expect normalized components, but this
  /// could be redundant work...
  void removeLowWeightComponents(std::vector<ComponentCache>& cmps) const {
    auto proj = [](auto& cmp) -> double& { return cmp.weight; };

    detail::normalizeWeights(cmps, proj);

    auto new_end = std::remove_if(cmps.begin(), cmps.end(), [&](auto& cmp) {
      return proj(cmp) < m_cfg.weightCutoff;
    });

    // In case we would remove all components, keep only the largest
    if (std::distance(cmps.begin(), new_end) == 0) {
      cmps = {*std::max_element(
          cmps.begin(), cmps.end(),
          [&](auto& a, auto& b) { return proj(a) < proj(b); })};
      cmps.front().weight = 1.0;
    } else {
      cmps.erase(new_end, cmps.end());
      detail::normalizeWeights(cmps, proj);
    }
  }

  /// Function that updates the stepper from the MultiTrajectory
  template <typename propagator_state_t, typename stepper_t>
  void updateStepper(propagator_state_t& state, const stepper_t& stepper,
                     const TemporaryStates& tmpStates) const {
    auto cmps = stepper.componentIterable(state.stepping);
    for (auto [idx, cmp] : zip(tmpStates.tips, cmps)) {
      // we set ignored components to missed, so we can remove them after
      // the loop
      if (tmpStates.weights.at(idx) < m_cfg.weightCutoff) {
        cmp.status() = IntersectionStatus::unreachable;
        continue;
      }

      auto proxy = tmpStates.traj.getTrackState(idx);

      cmp.pars() =
          MultiTrajectoryHelpers::freeFiltered(state.geoContext, proxy);
      cmp.cov() = proxy.filteredCovariance();
      cmp.weight() = tmpStates.weights.at(idx);
    }

    stepper.removeMissedComponents(state.stepping);

    // TODO we have two normalization passes here now, this can probably be
    // optimized
    detail::normalizeWeights(cmps,
                             [&](auto cmp) -> double& { return cmp.weight(); });
  }

  /// Function that updates the stepper from the ComponentCache
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  void updateStepper(propagator_state_t& state, const stepper_t& stepper,
                     const navigator_t& navigator,
                     const std::vector<ComponentCache>& componentCache) const {
    const auto& surface = *navigator.currentSurface(state.navigation);

    // Clear components before adding new ones
    stepper.clearComponents(state.stepping);

    // Finally loop over components
    for (const auto& [weight, pars, cov] : componentCache) {
      // Add the component to the stepper
      BoundTrackParameters bound(surface.getSharedPtr(), pars, cov,
                                 stepper.particleHypothesis(state.stepping));

      auto res = stepper.addComponent(state.stepping, std::move(bound), weight);

      if (!res.ok()) {
        ACTS_ERROR("Error adding component to MultiStepper");
        continue;
      }

      auto& cmp = *res;
      auto freeParams = cmp.pars();
      cmp.jacToGlobal() = surface.boundToFreeJacobian(
          state.geoContext, freeParams.template segment<3>(eFreePos0),
          freeParams.template segment<3>(eFreeDir0));
      cmp.pathAccumulated() = state.stepping.pathAccumulated;
      cmp.jacobian() = BoundMatrix::Identity();
      cmp.derivative() = FreeVector::Zero();
      cmp.jacTransport() = FreeMatrix::Identity();
    }
  }

  /// This function performs the kalman update, computes the new posterior
  /// weights, renormalizes all components, and does some statistics.
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  Result<void> kalmanUpdate(propagator_state_t& state, const stepper_t& stepper,
                            const navigator_t& navigator, result_type& result,
                            TemporaryStates& tmpStates,
                            const SourceLink& sourceLink) const {
    const auto& surface = *navigator.currentSurface(state.navigation);

    // Keep track of all created components for outlier handling
    std::vector<TrackIndexType> allTips;
    allTips.reserve(stepper.numberComponents(state.stepping));

    for (auto cmp : stepper.componentIterable(state.stepping)) {
      auto singleState = cmp.singleState(state);
      const auto& singleStepper = cmp.singleStepper(stepper);

      // Add a <mask> TrackState entry multi trajectory. This allocates storage
      // for all components, which we will set later.
      TrackStatePropMask mask =
          TrackStatePropMask::Predicted | TrackStatePropMask::Filtered |
          TrackStatePropMask::Jacobian | TrackStatePropMask::Calibrated;
      typename traj_t::TrackStateProxy trackStateProxy =
          tmpStates.traj.makeTrackState(mask, kTrackIndexInvalid);

      typename traj_t::ConstTrackStateProxy trackStateProxyConst{
          trackStateProxy};

      // Set the trackStateProxy components with the state from the ongoing
      // propagation
      {
        trackStateProxy.setReferenceSurface(surface.getSharedPtr());
        // Bind the transported state to the current surface
        auto res =
            singleStepper.boundState(singleState.stepping, surface, false);
        if (!res.ok()) {
          ACTS_DEBUG("Propagate to surface " << surface.geometryId()
                                             << " failed: " << res.error());
          return res.error();
        }
        const auto& [boundParams, jacobian, pathLength] = *res;

        // Fill the track state
        trackStateProxy.predicted() = boundParams.parameters();
        trackStateProxy.predictedCovariance() = singleState.stepping.cov;

        trackStateProxy.jacobian() = jacobian;
        trackStateProxy.pathLength() = pathLength;
      }

      // We have predicted parameters, so calibrate the uncalibrated input
      // measurement
      m_cfg.extensions.calibrator(state.geoContext, *m_cfg.calibrationContext,
                                  sourceLink, trackStateProxy);

      if (!m_cfg.extensions.outlierFinder(trackStateProxyConst)) {
        // Run Kalman update
        auto updateRes = m_cfg.extensions.updater(state.geoContext,
                                                  trackStateProxy, logger());
        if (!updateRes.ok()) {
          ACTS_DEBUG("Update step failed: " << updateRes.error());
          return updateRes.error();
        }

        tmpStates.tips.push_back(trackStateProxy.index());
        tmpStates.weights[trackStateProxy.index()] = cmp.weight();
      }

      allTips.push_back(trackStateProxy.index());
    }

    const bool isOutlier = tmpStates.tips.empty();

    if (!isOutlier) {
      computePosteriorWeights(tmpStates.traj, tmpStates.tips,
                              tmpStates.weights);
      detail::normalizeWeights(tmpStates.tips, [&](auto idx) -> double& {
        return tmpStates.weights.at(idx);
      });
    } else {
      auto cmps = stepper.componentIterable(state.stepping);
      for (const auto [cmp, idx] : zip(cmps, allTips)) {
        typename traj_t::TrackStateProxy trackStateProxy =
            tmpStates.traj.getTrackState(idx);

        // Set the filtered parameter index to be the same with predicted
        // parameter
        trackStateProxy.shareFrom(trackStateProxy,
                                  TrackStatePropMask::Predicted,
                                  TrackStatePropMask::Filtered);

        tmpStates.tips.push_back(trackStateProxy.index());
        tmpStates.weights[trackStateProxy.index()] = cmp.weight();
      }
    }

    // Do the statistics
    ++result.processedStates;
    if (!isOutlier) {
      ++result.measurementStates;
    }

    updateMultiTrajectory(result, tmpStates, surface,
                          isOutlier ? TrackStateType().setIsOutlier()
                                    : TrackStateType().setIsMeasurement());

    result.lastMeasurementTip = result.currentTip;
    result.lastMeasurementSurface = &surface;

    // Note, that we do not normalize the components here.
    // This must be done before initializing the backward pass.
    result.lastMeasurementComponents.clear();

    FiltProjector proj{tmpStates.traj, tmpStates.weights};
    for (const auto& idx : tmpStates.tips) {
      const auto& [w, p, c] = proj(idx);
      // TODO check why zero weight can occur
      if (w > 0.0) {
        result.lastMeasurementComponents.push_back({w, p, c});
      }
    }

    // Return success
    return Result<void>::success();
  }

  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  Result<void> noMeasurementUpdate(propagator_state_t& state,
                                   const stepper_t& stepper,
                                   const navigator_t& navigator,
                                   result_type& result,
                                   TemporaryStates& tmpStates,
                                   bool doCovTransport) const {
    const Surface& surface = *navigator.currentSurface(state.navigation);

    for (auto cmp : stepper.componentIterable(state.stepping)) {
      auto& singleState = cmp.state();
      const auto& singleStepper = cmp.singleStepper(stepper);

      // Add a <mask> TrackState entry multi trajectory. This allocates storage
      // for all components, which we will set later.
      TrackStatePropMask mask =
          TrackStatePropMask::Predicted | TrackStatePropMask::Jacobian;
      typename traj_t::TrackStateProxy trackStateProxy =
          tmpStates.traj.makeTrackState(mask, kTrackIndexInvalid);

      // Set the trackStateProxy components with the state from the ongoing
      // propagation
      {
        trackStateProxy.setReferenceSurface(surface.getSharedPtr());
        // Bind the transported state to the current surface
        auto res =
            singleStepper.boundState(singleState, surface, doCovTransport);
        if (!res.ok()) {
          return res.error();
        }
        const auto& [boundParams, jacobian, pathLength] = *res;

        // Fill the track state
        trackStateProxy.predicted() = boundParams.parameters();
        trackStateProxy.predictedCovariance() = singleState.cov;

        trackStateProxy.jacobian() = jacobian;
        trackStateProxy.pathLength() = pathLength;

        // Set the filtered parameter index to be the same with predicted
        // parameter
        trackStateProxy.shareFrom(trackStateProxy,
                                  TrackStatePropMask::Predicted,
                                  TrackStatePropMask::Filtered);
      }

      tmpStates.tips.push_back(trackStateProxy.index());
      tmpStates.weights[trackStateProxy.index()] = cmp.weight();
    }

    const bool precedingMeasurementExists = result.processedStates > 0;
    const bool isHole = surface.isSensitive();

    // Do the statistics
    ++result.processedStates;
    if (precedingMeasurementExists && isHole) {
      ++result.measurementHoles;
    }

    updateMultiTrajectory(
        result, tmpStates, surface,
        isHole ? TrackStateType().setIsHole() : TrackStateType());

    return Result<void>::success();
  }

  /// Apply the multiple scattering to the state
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  void applyMultipleScattering(propagator_state_t& state,
                               const stepper_t& stepper,
                               const navigator_t& navigator,
                               const MaterialUpdateMode& updateMode) const {
    for (auto cmp : stepper.componentIterable(state.stepping)) {
      auto singleState = cmp.singleState(state);
      const auto& singleStepper = cmp.singleStepper(stepper);

      detail::performMaterialInteraction(
          singleState, singleStepper, navigator, updateMode,
          NoiseUpdateMode::addNoise, m_cfg.multipleScattering, false, logger());

      assert(singleState.stepping.cov.array().isFinite().all() &&
             "covariance not finite after multi scattering");
    }
  }

  void updateMultiTrajectory(result_type& result,
                             const TemporaryStates& tmpStates,
                             const Surface& surface,
                             TrackStateType type) const {
    using PrtProjector =
        MultiTrajectoryProjector<StatesType::ePredicted, traj_t>;
    using FltProjector =
        MultiTrajectoryProjector<StatesType::eFiltered, traj_t>;

    if (!m_cfg.inReversePass) {
      const bool isMeasurement = type.isMeasurement();
      const bool isOutlier = type.isOutlier();
      const bool isHole = type.isHole();

      const bool hasMaterial = surface.surfaceMaterial() != nullptr;
      const auto firstCmpProxy =
          tmpStates.traj.getTrackState(tmpStates.tips.front());

      auto mask = TrackStatePropMask::Predicted | TrackStatePropMask::Smoothed;
      if (isMeasurement) {
        mask |= TrackStatePropMask::Calibrated | TrackStatePropMask::Filtered;
      } else if (isOutlier) {
        mask |= TrackStatePropMask::Calibrated;
      }

      auto proxy = result.fittedStates->makeTrackState(mask, result.currentTip);
      result.currentTip = proxy.index();

      proxy.setReferenceSurface(surface.getSharedPtr());
      if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Calibrated)) {
        // copy source link, calibrated measurement, and subspace
        proxy.copyFrom(firstCmpProxy, TrackStatePropMask::Calibrated);
      }

      auto [prtMean, prtCov] =
          mergeGaussianMixture(tmpStates.tips, surface, m_cfg.mergeMethod,
                               PrtProjector{tmpStates.traj, tmpStates.weights});
      proxy.predicted() = prtMean;
      proxy.predictedCovariance() = prtCov;

      proxy.typeFlags().setHasParameters();
      proxy.typeFlags().setHasMaterial(hasMaterial);

      if (isMeasurement) {
        auto [fltMean, fltCov] = mergeGaussianMixture(
            tmpStates.tips, surface, m_cfg.mergeMethod,
            FltProjector{tmpStates.traj, tmpStates.weights});
        proxy.filtered() = fltMean;
        proxy.filteredCovariance() = fltCov;

        proxy.typeFlags().setIsMeasurement();
      } else {
        proxy.shareFrom(TrackStatePropMask::Predicted,
                        TrackStatePropMask::Filtered);
      }

      proxy.typeFlags().setIsOutlier(isOutlier);
      proxy.typeFlags().setIsHole(isHole);

      proxy.smoothed() = BoundVector::Constant(-2);
      proxy.smoothedCovariance() = BoundSquareMatrix::Constant(-2);

    } else {
      assert((result.currentTip != kTrackIndexInvalid && "tip not valid"));

      result.fittedStates->applyBackwards(
          result.currentTip, [&](auto trackState) {
            auto fSurface = &trackState.referenceSurface();
            if (fSurface == &surface) {
              result.surfacesVisitedBwdAgain.push_back(&surface);

              if (trackState.hasSmoothed()) {
                const auto [smtMean, smtCov] = mergeGaussianMixture(
                    tmpStates.tips, surface, m_cfg.mergeMethod,
                    FltProjector{tmpStates.traj, tmpStates.weights});

                trackState.smoothed() = smtMean;
                trackState.smoothedCovariance() = smtCov;
              }
              return false;
            }
            return true;
          });
    }
  }

  /// Set the relevant options that can be set from the Options struct all in
  /// one place
  void setOptions(const GsfOptions<traj_t>& options) {
    m_cfg.maxComponents = options.maxComponents;
    m_cfg.extensions = options.extensions;
    m_cfg.abortOnError = options.abortOnError;
    m_cfg.disableAllMaterialHandling = options.disableAllMaterialHandling;
    m_cfg.weightCutoff = options.weightCutoff;
    m_cfg.mergeMethod = options.componentMergeMethod;
    m_cfg.calibrationContext = &options.calibrationContext.get();
  }
};

}  // namespace Acts::detail

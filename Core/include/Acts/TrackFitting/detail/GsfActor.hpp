// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiComponentBoundTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GsfError.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/BetheHeitlerApprox.hpp"
#include "Acts/TrackFitting/detail/GsfSmoothing.hpp"
#include "Acts/TrackFitting/detail/GsfUtils.hpp"
#include "Acts/TrackFitting/detail/KLMixtureReduction.hpp"
#include "Acts/TrackFitting/detail/KalmanUpdateHelpers.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <ios>
#include <map>
#include <numeric>

namespace Acts {
namespace detail {

template <typename traj_t>
struct GsfResult {
  /// The multi-trajectory which stores the graph of components
  std::shared_ptr<traj_t> fittedStates;

  /// This provides the weights for the states in the MultiTrajectory. Each
  /// entry maps to one track state. TODO This is a workaround until the
  /// MultiTrajectory can handle weights
  std::map<MultiTrajectoryTraits::IndexType, ActsScalar> weightsOfStates;

  /// The current indexes for the newest components in the multi trajectory
  /// (this includes material, hole and outlier states)
  std::vector<MultiTrajectoryTraits::IndexType> currentTips;

  /// The last tips referring to a measuerement state so we do not need so
  /// search them recursively later
  std::vector<MultiTrajectoryTraits::IndexType> lastMeasurementTips;

  /// We must capture the parent tips to ensure that we can keep track of the
  /// last states in the multitrajectory after the component convolution and
  /// reduction
  std::vector<MultiTrajectoryTraits::IndexType> parentTips;

  /// Some counting
  std::size_t measurementStates = 0;
  std::size_t measurementHoles = 0;
  std::size_t processedStates = 0;
  std::set<Acts::GeometryIdentifier> visitedSurfaces;
  std::vector<const Acts::Surface*> missedActiveSurfaces;

  // Propagate potential errors to the outside
  Result<void> result{Result<void>::success()};

  // Used for workaround to initialize MT correctly
  bool haveInitializedResult = false;
};

/// The actor carrying out the GSF algorithm
template <typename bethe_heitler_approx_t, typename traj_t>
struct GsfActor {
  /// Enforce default construction
  GsfActor() = default;

  /// Broadcast the result_type
  using result_type = GsfResult<traj_t>;

  // Actor configuration
  struct Config {
    /// Maximum number of components which the GSF should handle
    std::size_t maxComponents = 16;

    /// Input measurements
    std::map<GeometryIdentifier, std::reference_wrapper<const SourceLink>>
        inputMeasurements;

    /// Bethe Heitler Approximator pointer. The fitter holds the approximator
    /// instance TODO if we somehow could initialize a reference here...
    const bethe_heitler_approx_t* bethe_heitler_approx = nullptr;

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// When to discard components
    double weightCutoff = 1.0e-4;

    /// A not so nice workaround to get the first backward state in the
    /// MultiTrajectory for the DirectNavigator
    std::function<void(result_type&, const LoggerWrapper&)> resultInitializer;

    /// When this option is enabled, material information on all surfaces is
    /// ignored. This disables the component convolution as well as the handling
    /// of energy. This may be useful for debugging.
    bool disableAllMaterialHandling = false;

    /// Whether to abort immediately when an error occurs
    bool abortOnError = false;

    /// The extensions
    GsfExtensions<traj_t> extensions;
  } m_cfg;

  /// Stores meta information about the components
  struct MetaCache {
    /// Where to find the parent component in the MultiTrajectory
    MultiTrajectoryTraits::IndexType parentIndex;

    /// Other quantities TODO are they really needed here? seems they are
    /// reinitialized to Identity etc.
    BoundMatrix jacobian;
    BoundToFreeMatrix jacToGlobal;
    FreeMatrix jacTransport;
    FreeVector derivative;

    /// We need to preserve the path length
    ActsScalar pathLength;
  };

  /// Stores parameters of a gaussian component
  struct ParameterCache {
    ActsScalar weight;
    BoundVector boundPars;
    std::optional<BoundSymMatrix> boundCov;
  };

  /// Broadcast Cache Type
  using ComponentCache = std::tuple<ParameterCache, MetaCache>;

  /// @brief GSF actor operation
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
    assert(result.fittedStates && "No MultiTrajectory set");
    const auto& logger = state.options.logger;

    // Prints some VERBOSE things and performs some asserts. Can be removed
    // without change of behaviour
    const detail::ScopedGsfInfoPrinterAndChecker printer(state, stepper);

    // Set error or abort utility
    auto set_error_or_abort = [&](auto error) {
      if (m_cfg.abortOnError) {
        std::abort();
      } else {
        result.result = error;
      }
    };

    // Count the states of the components, this is necessary to evaluate if
    // really all components are on a surface TODO Not sure why this is not
    // garantueed by having currentSurface pointer set
    const auto [missed_count, reachable_count] = [&]() {
      std::size_t missed = 0;
      std::size_t reachable = 0;
      for (auto cmp : stepper.constComponentIterable(state.stepping)) {
        using Status = Acts::Intersection3D::Status;

        // clang-format off
          switch (cmp.status()) {
            break; case Status::missed: ++missed;
            break; case Status::reachable: ++reachable;
            break; default: {}
          }
        // clang-format on
      }
      return std::make_tuple(missed, reachable);
    }();

    // Workaround to initialize e.g. MultiTrajectory in backward mode
    if (!result.haveInitializedResult && m_cfg.resultInitializer) {
      m_cfg.resultInitializer(result, logger);
      result.haveInitializedResult = true;
    }

    // Initialize the tips if they are empty (should only happen at first pass)
    if (result.parentTips.empty()) {
      result.parentTips.resize(stepper.numberComponents(state.stepping),
                               MultiTrajectoryTraits::kInvalid);
    }

    if (result.parentTips.size() != stepper.numberComponents(state.stepping)) {
      ACTS_ERROR("component number mismatch:"
                 << result.parentTips.size() << " vs "
                 << stepper.numberComponents(state.stepping));

      return set_error_or_abort(GsfError::ComponentNumberMismatch);
    }

    // There seem to be cases where this is not always after initializing the
    // navigation from a surface. Some later functions assume this criterium
    // to be fulfilled.
    bool on_surface = reachable_count == 0 &&
                      missed_count < stepper.numberComponents(state.stepping);

    // We only need to do something if we are on a surface
    if (state.navigation.currentSurface && on_surface) {
      const auto& surface = *state.navigation.currentSurface;
      ACTS_VERBOSE("Step is at surface " << surface.geometryId());

      // Early return if we already were on this surface TODO why is this
      // necessary
      const auto [it, success] =
          result.visitedSurfaces.insert(surface.geometryId());

      if (!success) {
        ACTS_VERBOSE("Already visited surface, return");
        return;
      }

      removeMissedComponents(state, stepper, result.parentTips);

      // Check what we have on this surface
      const auto found_source_link =
          m_cfg.inputMeasurements.find(surface.geometryId());
      const bool haveMaterial =
          state.navigation.currentSurface->surfaceMaterial() &&
          !m_cfg.disableAllMaterialHandling;
      const bool haveMeasurement =
          found_source_link != m_cfg.inputMeasurements.end();

      ACTS_VERBOSE(std::boolalpha << "haveMaterial " << haveMaterial
                                  << ", haveMeasurement: " << haveMeasurement);

      ////////////////////////
      // The Core Algorithm
      ////////////////////////

      // Early return if nothing happens
      if (not haveMaterial && not haveMeasurement) {
        // No hole before first measurement
        if (result.processedStates > 0 && surface.associatedDetectorElement()) {
          noMeasurementUpdate(state, stepper, result, true);
          result.parentTips = result.currentTips;
        }
        return;
      }

      for (auto cmp : stepper.componentIterable(state.stepping)) {
        auto singleState = cmp.singleState(state);
        cmp.singleStepper(stepper).transportCovarianceToBound(
            singleState.stepping, surface);
      }

      if (haveMaterial) {
        if (haveMeasurement) {
          applyMultipleScattering(state, stepper,
                                  MaterialUpdateStage::PreUpdate);
        } else {
          applyMultipleScattering(state, stepper,
                                  MaterialUpdateStage::FullUpdate);
        }
      }

      // We do not need the component cache here, we can just update our stepper
      // state with the filtered components.
      // NOTE because of early return before we know that we have a measurement
      if (not haveMaterial) {
        kalmanUpdate(state, stepper, result, found_source_link->second);

        result.parentTips = updateStepper(state, stepper, result);

      }
      // We have material, we thus need a component cache since we will
      // convolute the components and later reduce them again before updating
      // the stepper
      else {
        std::vector<ComponentCache> componentCache;

        if (haveMeasurement) {
          kalmanUpdate(state, stepper, result, found_source_link->second);

          convoluteComponents(state, stepper, result, componentCache);
        } else {
          noMeasurementUpdate(state, stepper, result, false);

          convoluteComponents(state, stepper, result, componentCache);
        }

        reduceComponents(stepper, surface, componentCache);

        removeLowWeightComponents(componentCache);

        result.parentTips = updateStepper(state, stepper, componentCache);
      }

      // If we only done preUpdate before, now do postUpdate
      if (haveMaterial && haveMeasurement) {
        applyMultipleScattering(state, stepper,
                                MaterialUpdateStage::PostUpdate);
      }
    }
  }

  template <typename propagator_state_t, typename stepper_t>
  void convoluteComponents(propagator_state_t& state, const stepper_t& stepper,
                           const result_type& result,
                           std::vector<ComponentCache>& componentCache) const {
    auto cmps = stepper.componentIterable(state.stepping);
    for (auto [idx, cmp] : zip(result.currentTips, cmps)) {
      auto proxy = result.fittedStates->getTrackState(idx);

      MetaCache mcache;
      mcache.parentIndex = idx;
      mcache.jacobian = cmp.jacobian();
      mcache.jacToGlobal = cmp.jacToGlobal();
      mcache.jacTransport = cmp.jacTransport();
      mcache.derivative = cmp.derivative();
      mcache.pathLength = cmp.pathAccumulated();

      BoundTrackParameters bound(proxy.referenceSurface().getSharedPtr(),
                                 proxy.filtered(), proxy.filteredCovariance());

      applyBetheHeitler(state, bound, result.weightsOfStates.at(idx), mcache,
                        componentCache);
    }
  }

  template <typename propagator_state_t>
  void applyBetheHeitler(const propagator_state_t& state,
                         const BoundTrackParameters& old_bound,
                         const double old_weight, const MetaCache& metaCache,
                         std::vector<ComponentCache>& componentCaches) const {
    const auto& logger = state.options.logger;
    const auto& surface = *state.navigation.currentSurface;
    const auto p_prev = old_bound.absoluteMomentum();

    // Evaluate material slab
    auto slab = surface.surfaceMaterial()->materialSlab(
        old_bound.position(state.stepping.geoContext), state.stepping.navDir,
        MaterialUpdateStage::FullUpdate);

    auto pathCorrection =
        surface.pathCorrection(state.stepping.geoContext,
                               old_bound.position(state.stepping.geoContext),
                               old_bound.unitDirection());
    slab.scaleThickness(pathCorrection);

    // Get the mixture
    const auto mixture =
        m_cfg.bethe_heitler_approx->mixture(slab.thicknessInX0());

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
        if (state.stepping.navDir == NavigationDirection::Forward)
          return p_prev * (gaussian.mean - 1.);
        else
          return p_prev * (1. / gaussian.mean - 1.);
      }();

      throw_assert(p_prev + delta_p > 0.,
                   "new momentum after bethe-heitler must be > 0, p_prev= "
                       << p_prev << ", delta_p=" << delta_p
                       << ", gaussian mean: " << gaussian.mean);
      new_pars[eBoundQOverP] = old_bound.charge() / (p_prev + delta_p);

      // compute inverse variance of p from mixture and update covariance
      auto new_cov = std::move(old_bound.covariance());

      if (new_cov.has_value()) {
        const auto varInvP = [&]() {
          if (state.stepping.navDir == NavigationDirection::Forward) {
            const auto f = 1. / (p_prev * gaussian.mean);
            return f * f * gaussian.var;
          } else {
            return gaussian.var / (p_prev * p_prev);
          }
        }();

        (*new_cov)(eBoundQOverP, eBoundQOverP) += varInvP;
        throw_assert(
            std::isfinite((*new_cov)(eBoundQOverP, eBoundQOverP)),
            "cov not finite, varInvP=" << varInvP << ", p_prev=" << p_prev
                                       << ", gaussian.mean=" << gaussian.mean
                                       << ", gaussian.var=" << gaussian.var);
      }

      // Set the remaining things and push to vector
      componentCaches.push_back(
          {ParameterCache{new_weight, new_pars, new_cov}, metaCache});
    }
  }

  template <typename stepper_t>
  void reduceComponents(const stepper_t& stepper, const Surface& surface,
                        std::vector<ComponentCache>& cmps) const {
    // Final component number
    const auto final_cmp_number = std::min(
        static_cast<std::size_t>(stepper.maxComponents), m_cfg.maxComponents);

    auto proj = [](auto& a) -> decltype(auto) { return std::get<0>(a); };

    // We must differ between surface types, since there can be different
    // local coordinates
    // TODO add other surface types
    switch (surface.type()) {
      case Surface::Cylinder: {
        // The cylinder coordinate is phi*R, so we need to pass R
        detail::AngleDescription::Cylinder angle_desc;
        std::get<0>(angle_desc).constant =
            static_cast<const CylinderSurface&>(surface).bounds().get(
                CylinderBounds::eR);

        detail::reduceWithKLDistance(cmps, final_cmp_number, proj, angle_desc);
      } break;
      default: {
        detail::reduceWithKLDistance(cmps, final_cmp_number, proj);
      }
    }
  }

  /// Removes the components which are missed and update the list of parent tips
  /// for the MultiTrajectory
  template <typename propagator_state_t, typename stepper_t>
  void removeMissedComponents(
      propagator_state_t& state, const stepper_t& stepper,
      std::vector<MultiTrajectoryTraits::IndexType>& tips) const {
    std::vector<MultiTrajectoryTraits::IndexType> new_tips;
    auto components = stepper.componentIterable(state.stepping);
    double sum_w = 0.0;

    for (auto [tip, cmp] : zip(tips, components)) {
      if (cmp.status() == Intersection3D::Status::onSurface) {
        sum_w += cmp.weight();
        new_tips.push_back(tip);
      }
    }

    // If the remaining weights are close to zero, re-sanitize all weights
    if (sum_w < m_cfg.weightCutoff) {
      for (auto cmp : components) {
        cmp.weight() = 1.0;
      }
    }

    tips = new_tips;
    stepper.removeMissedComponents(state.stepping);

    detail::normalizeWeights(components,
                             [](auto& cmp) -> double& { return cmp.weight(); });

    throw_assert(stepper.numberComponents(state.stepping) == tips.size(),
                 "size mismatch");
  }

  /// Remove components with low weights and renormalize from the component
  /// cache
  /// TODO This function does not expect normalized components, but this
  /// could be redundant work...
  void removeLowWeightComponents(std::vector<ComponentCache>& cmps) const {
    auto proj = [](auto& cmp) -> double& { return std::get<0>(cmp).weight; };

    detail::normalizeWeights(cmps, proj);

    auto new_end = std::remove_if(cmps.begin(), cmps.end(), [&](auto& cmp) {
      return proj(cmp) < m_cfg.weightCutoff;
    });
    cmps.erase(new_end, cmps.end());

    detail::normalizeWeights(cmps, proj);
  }

  /// Function that updates the stepper from the MultiTrajectory
  template <typename propagator_state_t, typename stepper_t>
  std::vector<MultiTrajectoryTraits::IndexType> updateStepper(
      propagator_state_t& state, const stepper_t& stepper,
      const result_type& result) const {
    std::vector<MultiTrajectoryTraits::IndexType> newTips;

    auto cmps = stepper.componentIterable(state.stepping);

    for (auto [idx, cmp] : zip(result.currentTips, cmps)) {
      // we set ignored components to missed, so we can remove them after
      // the loop
      if (result.weightsOfStates.at(idx) < m_cfg.weightCutoff) {
        cmp.status() = Intersection3D::Status::missed;
        continue;
      }

      auto proxy = result.fittedStates->getTrackState(idx);

      cmp.pars() =
          MultiTrajectoryHelpers::freeFiltered(state.options.geoContext, proxy);
      cmp.cov() = proxy.filteredCovariance();
      cmp.weight() = result.weightsOfStates.at(idx);

      newTips.push_back(idx);
    }

    stepper.removeMissedComponents(state.stepping);

    // TODO we have two normalization passes here now, this can probably be
    // optimized
    detail::normalizeWeights(cmps,
                             [&](auto cmp) -> double& { return cmp.weight(); });

    return newTips;
  }

  /// Function that updates the stepper from the ComponentCache
  template <typename propagator_state_t, typename stepper_t>
  std::vector<MultiTrajectoryTraits::IndexType> updateStepper(
      propagator_state_t& state, const stepper_t& stepper,
      const std::vector<ComponentCache>& componentCache) const {
    const auto& surface = *state.navigation.currentSurface;
    const auto& logger = state.options.logger;

    // We collect new tips in the loop
    std::vector<MultiTrajectoryTraits::IndexType> new_parent_tips;

    // Clear components before adding new ones
    stepper.clearComponents(state.stepping);

    // Finally loop over components
    for (const auto& [pcache, meta] : componentCache) {
      const auto& [weight, pars, cov] = pcache;

      // Add the component to the stepper
      const BoundTrackParameters bound(surface.getSharedPtr(), pars, cov);

      auto res = stepper.addComponent(state.stepping, std::move(bound), weight);

      if (!res.ok()) {
        ACTS_ERROR("Error adding component to MultiStepper");
        continue;
      }

      // Only add index after we are sure the component was successful added
      new_parent_tips.push_back(meta.parentIndex);

      auto& cmp = *res;
      cmp.jacobian() = meta.jacobian;
      cmp.jacToGlobal() = meta.jacToGlobal;
      cmp.pathAccumulated() = meta.pathLength;
      cmp.derivative() = meta.derivative;
      cmp.jacTransport() = meta.jacTransport;
    }

    return new_parent_tips;
  }

  /// This function performs the kalman update, computes the new posterior
  /// weights, renormalizes all components, and does some statistics.
  template <typename propagator_state_t, typename stepper_t>
  Result<void> kalmanUpdate(propagator_state_t& state, const stepper_t& stepper,
                            result_type& result,
                            const SourceLink& source_link) const {
    const auto& surface = *state.navigation.currentSurface;

    // We will overwrite this soon with new components
    result.currentTips.clear();

    // Boolean flag, to distinguish measurement and outlier states. This flag
    // is only modified by the valid-measurement-branch, so only if there
    // isn't any valid measuerement state, the flag stays false and the state
    // is thus counted as an outlier
    bool is_valid_measurement = false;

    auto cmps = stepper.componentIterable(state.stepping);
    for (auto [idx, cmp] : zip(result.parentTips, cmps)) {
      auto singleState = cmp.singleState(state);
      const auto& singleStepper = cmp.singleStepper(stepper);

      auto trackStateProxyRes = detail::kalmanHandleMeasurement(
          singleState, singleStepper, m_cfg.extensions, surface, source_link,
          *result.fittedStates, idx, false);

      if (!trackStateProxyRes.ok()) {
        return trackStateProxyRes.error();
      }

      const auto& trackStateProxy = *trackStateProxyRes;

      // If at least one component is no outlier, we consider the whole thing
      // as a measuerementState
      if (trackStateProxy.typeFlags().test(
              Acts::TrackStateFlag::MeasurementFlag)) {
        is_valid_measurement = true;
      }

      result.currentTips.push_back(trackStateProxy.index());
      result.weightsOfStates[result.currentTips.back()] = cmp.weight();
    }

    computePosteriorWeights(*result.fittedStates, result.currentTips,
                            result.weightsOfStates);

    detail::normalizeWeights(result.currentTips, [&](auto idx) -> double& {
      return result.weightsOfStates.at(idx);
    });

    // Do the statistics
    ++result.processedStates;

    // We also need to save outlier states here, otherwise they would not be
    // included in the MT if they are at the end of the track
    result.lastMeasurementTips = result.currentTips;

    // TODO should outlier states also be counted here?
    if (is_valid_measurement) {
      ++result.measurementStates;
    }

    // Return sucess
    return Acts::Result<void>::success();
  }

  template <typename propagator_state_t, typename stepper_t>
  Result<void> noMeasurementUpdate(propagator_state_t& state,
                                   const stepper_t& stepper,
                                   result_type& result,
                                   bool doCovTransport) const {
    const auto& surface = *state.navigation.currentSurface;

    // We will overwrite this soon with new components
    result.currentTips.clear();

    // Initialize as true, so that any component can flip it. However, all
    // components should behave the same
    bool is_hole = true;

    auto cmps = stepper.componentIterable(state.stepping);
    for (auto [idx, cmp] : zip(result.parentTips, cmps)) {
      auto singleState = cmp.singleState(state);
      const auto& singleStepper = cmp.singleStepper(stepper);

      // There is some redundant checking inside this function, but do this for
      // now until we measure this is significant
      auto trackStateProxyRes = detail::kalmanHandleNoMeasurement(
          singleState, singleStepper, surface, *result.fittedStates, idx,
          doCovTransport);

      if (!trackStateProxyRes.ok()) {
        return trackStateProxyRes.error();
      }

      const auto& trackStateProxy = *trackStateProxyRes;

      if (not trackStateProxy.typeFlags().test(TrackStateFlag::HoleFlag)) {
        is_hole = false;
      }

      result.currentTips.push_back(trackStateProxy.index());
      result.weightsOfStates[result.currentTips.back()] = cmp.weight();
    }

    // These things should only be done once for all components
    if (is_hole) {
      result.missedActiveSurfaces.push_back(&surface);
      ++result.measurementHoles;
    }

    ++result.processedStates;

    return Result<void>::success();
  }

  /// Apply the multipe scattering to the state
  template <typename propagator_state_t, typename stepper_t>
  void applyMultipleScattering(propagator_state_t& state,
                               const stepper_t& stepper,
                               const MaterialUpdateStage& updateStage =
                                   MaterialUpdateStage::FullUpdate) const {
    const auto& logger = state.options.logger;
    const auto& surface = *state.navigation.currentSurface;

    for (auto cmp : stepper.componentIterable(state.stepping)) {
      auto singleState = cmp.singleState(state);
      const auto& singleStepper = cmp.singleStepper(stepper);

      detail::PointwiseMaterialInteraction interaction(&surface, singleState,
                                                       singleStepper);
      if (interaction.evaluateMaterialSlab(singleState, updateStage)) {
        // In the Gsf we only need to handle the multiple scattering
        interaction.evaluatePointwiseMaterialInteraction(
            m_cfg.multipleScattering, false);

        // Screen out material effects info
        ACTS_VERBOSE("Material effects on surface: "
                     << surface.geometryId()
                     << " at update stage: " << updateStage << " are :");
        ACTS_VERBOSE("eLoss = "
                     << interaction.Eloss << ", "
                     << "variancePhi = " << interaction.variancePhi << ", "
                     << "varianceTheta = " << interaction.varianceTheta << ", "
                     << "varianceQoverP = " << interaction.varianceQoverP);

        // Update the state and stepper with material effects
        interaction.updateState(singleState, singleStepper, addNoise);

        throw_assert(singleState.stepping.cov.array().isFinite().all(),
                     "covariance not finite after update");
      }
    }
  }
};

}  // namespace detail
}  // namespace Acts

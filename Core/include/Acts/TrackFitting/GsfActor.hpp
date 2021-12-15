// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
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
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/BetheHeitlerApprox.hpp"
#include "Acts/TrackFitting/GsfError.hpp"
#include "Acts/TrackFitting/detail/GsfSmoothing.hpp"
#include "Acts/TrackFitting/detail/GsfUtils.hpp"
#include "Acts/TrackFitting/detail/KLMixtureReduction.hpp"
#include "Acts/Utilities/Overload.hpp"

#include <ios>
#include <map>
#include <numeric>


#define RETURN_ERROR_OR_ABORT_ACTOR(error) \
  if (m_cfg.abortOnError) {                \
    std::abort();                          \
  } else {                                 \
    return error;                          \
  }

#define SET_ERROR_AND_RETURN_OR_ABORT_ACTOR(error) \
  if (m_cfg.abortOnError) {                        \
    std::abort();                                  \
  } else {                                         \
    result.result = error;                         \
    return;                                        \
  }

namespace Acts {

struct GsfResult {
  /// The multi-trajectory which stores the graph of components
  MultiTrajectory fittedStates;
  std::map<size_t, ActsScalar> weightsOfStates;

  /// The current indexes for the newest components in the multi trajectory
  std::vector<std::size_t> currentTips;

  /// The indices of the parent states in the multi trajectory. This is needed
  /// since there can be cases with splitting but no Kalman-Update, so we need
  /// to preserve this information
  std::vector<std::size_t> parentTips;

  /// The number of measurement states created
  std::size_t measurementStates = 0;
  std::size_t processedStates = 0;
  std::set<Acts::GeometryIdentifier> visitedSurfaces;

  // Propagate potential errors to the outside
  Result<void> result{Result<void>::success()};

  // Used for workaround to initialize MT correctly
  bool haveInitializedMT = false;
};

/// The actor carrying out the GSF algorithm
struct GsfActor {
  /// Enforce default construction
  GsfActor() = default;

  /// Broadcast the result_type
  using result_type = GsfResult;

  // Actor configuration
  struct Config {
    /// Maximum number of components which the GSF should handle
    std::size_t maxComponents = 16;

    /// Input measurements
    std::map<GeometryIdentifier, std::reference_wrapper<const SourceLink>>
        inputMeasurements;

    /// Bethe Heitler Approximator pointer. The fitter holds the approximator
    /// instance TODO if we somehow could initialize a reference here...
    const detail::BHApprox* bethe_heitler_approx = nullptr;

    /// Number of processed states
    std::size_t processedStates = 0;

    /// Allow to configure surfaces to skip
    std::set<GeometryIdentifier> surfacesToSkip;

    /// Wether to transport covariance
    bool doCovTransport = true;

    /// Whether to consider multiple scattering.
    bool multipleScattering = true;

    /// Whether to consider energy loss.
    bool energyLoss = true;

    /// Whether to abort immediately when an error occurs
    bool abortOnError = false;

    /// When to discard components
    double weightCutoff = 1.0e-4;

    /// A not so nice workaround to get the first backward state in the
    /// MultiTrajectory for the DirectNavigator
    std::function<void(result_type&, const LoggerWrapper&)>
        multiTrajectoryInitializer;

    /// We can disable component splitting for debugging or so
    bool applyMaterialEffects = true;
    
    /// The extensions
    KalmanFitterExtensions extensions;
  } m_cfg;

  /// Broadcast Cache Type
  using TrackProxy = typename MultiTrajectory::TrackStateProxy;
  using ComponentCache =
      std::tuple<std::variant<detail::GsfComponentParameterCache, TrackProxy>,
                 detail::GsfComponentMetaCache>;

  struct ParametersCacheProjector {
    auto& operator()(ComponentCache& cache) const {
      return std::get<detail::GsfComponentParameterCache>(std::get<0>(cache));
    }
    const auto& operator()(const ComponentCache& cache) const {
      return std::get<detail::GsfComponentParameterCache>(std::get<0>(cache));
    }
  };

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
    const auto& logger = state.options.logger;

    throw_assert(detail::weightsAreNormalized(
                     stepper.constComponentIterable(state.stepping),
                     [](const auto& cmp) { return cmp.weight(); }),
                 "not normalized at start of operator()");

    // A class that prints information about the state on construction and
    // destruction
    class ScopedInfoPrinter {
      const propagator_state_t& m_state;
      const stepper_t& m_stepper;
      double m_p_initial;

      const auto& logger() const { return m_state.options.logger(); }

      void print_component_stats() const {
        std::size_t i = 0;
        for (auto cmp : m_stepper.constComponentIterable(m_state.stepping)) {
          auto getVector = [&](auto idx) {
            return cmp.pars().template segment<3>(idx).transpose();
          };
          ACTS_VERBOSE("  #" << i++ << " pos: " << getVector(eFreePos0)
                             << ", dir: " << getVector(eFreeDir0)
                             << ", weight: " << cmp.weight()
                             << ", status: " << cmp.status()
                             << ", qop: " << cmp.pars()[eFreeQOverP]);
        }
      }

     public:
      ScopedInfoPrinter(const propagator_state_t& state,
                        const stepper_t& stepper)
          : m_state(state),
            m_stepper(stepper),
            m_p_initial(stepper.momentum(state.stepping)) {
        // Some initial printing
        ACTS_VERBOSE("Gsf step "
                     << state.stepping.steps << " at mean position "
                     << stepper.position(state.stepping).transpose()
                     << " with direction "
                     << stepper.direction(state.stepping).transpose()
                     << " and momentum " << stepper.momentum(state.stepping)
                     << " and charge " << stepper.charge(state.stepping));
        ACTS_VERBOSE(
            "Propagation is in "
            << (state.stepping.navDir == forward ? "forward" : "backward")
            << " mode");
        print_component_stats();
      }

      ~ScopedInfoPrinter() {
        if (m_state.navigation.currentSurface) {
          const auto p_final = m_stepper.momentum(m_state.stepping);
          ACTS_VERBOSE("Component status at end of step:");
          print_component_stats();
          ACTS_VERBOSE("Delta Momentum = " << std::setprecision(5)
                                           << p_final - m_p_initial);
        }
      }
    };

    const ScopedInfoPrinter printer(state, stepper);

    // Count the states of the components, this is necessary to evaluate if
    // really all components are on a surface TODO Not sure why this is not
    // garantueed by having currentSurface pointer set
    const auto [missed_count, reachable_count] = [&]() {
      std::size_t missed = 0;
      std::size_t reachable = 0;
      for (auto cmp : stepper.componentIterable(state.stepping)) {
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

    // Workaround to initialize MT in backward mode
    if (!result.haveInitializedMT && m_cfg.multiTrajectoryInitializer) {
      m_cfg.multiTrajectoryInitializer(result, logger);
      result.haveInitializedMT = true;
    }

    // Initialize current tips on first pass
    if (result.parentTips.empty()) {
      result.parentTips.resize(stepper.numberComponents(state.stepping),
                               SIZE_MAX);
    }

    // Check if we have the right number of components
    if (result.parentTips.size() != stepper.numberComponents(state.stepping)) {
      ACTS_ERROR("component number mismatch:"
                 << result.parentTips.size() << " vs "
                 << stepper.numberComponents(state.stepping));

      SET_ERROR_AND_RETURN_OR_ABORT_ACTOR(GsfError::ComponentNumberMismatch);
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

      // Early return if the surfaces should be skipped
      if (m_cfg.surfacesToSkip.find(surface.geometryId()) !=
          m_cfg.surfacesToSkip.end()) {
        ACTS_VERBOSE("Surface is configured to be skipped");
        return;
      }

      // Remove missed surfaces and adjust momenta
      removeMissedComponents(state, stepper, result.parentTips);
      throw_assert(
          result.parentTips.size() == stepper.numberComponents(state.stepping),
          "size mismatch (parentTips="
              << result.parentTips.size()
              << ", nCmps=" << stepper.numberComponents(state.stepping));

      // Check what we have on this surface
      const auto found_source_link =
          m_cfg.inputMeasurements.find(surface.geometryId());
      const bool haveMaterial =
          state.navigation.currentSurface->surfaceMaterial() &&
          m_cfg.applyMaterialEffects;
      const bool haveMeasurement =
          found_source_link != m_cfg.inputMeasurements.end();

      ACTS_VERBOSE(std::boolalpha << "haveMaterial " << haveMaterial
                                  << ", haveMeasurement: " << haveMeasurement);

      // Early return if nothing happens
      if (not haveMaterial && not haveMeasurement) {
        return;
      }

      // Create Cache
      thread_local std::vector<ComponentCache> componentCache;
      componentCache.clear();

      // Projectors
      auto mapToProxyAndWeight = [&](auto& cmp) {
        auto& proxy = std::get<TrackProxy>(std::get<0>(cmp));
        return std::tie(proxy, result.weightsOfStates.at(proxy.index()));
      };

      auto mapParsCacheToWeightParsCov = [&](auto& variant) -> decltype(auto) {
        return std::get<detail::GsfComponentParameterCache>(variant);
      };

      auto mapProxyToWeightParsCov = [&](auto& variant) {
        auto& proxy = std::get<TrackProxy>(variant);
        return std::make_tuple(result.weightsOfStates.at(proxy.index()),
                               proxy.filtered(), proxy.filteredCovariance());
      };

      auto mapToWeight = [&](auto& cmp) -> decltype(auto) {
        constexpr bool C =
            std::is_const_v<std::remove_reference_t<decltype(cmp)>>;

        using T1 = detail::GsfComponentParameterCache;
        using T2 = MultiTrajectory::TrackStateProxy;
        using R = std::conditional_t<C, double, double&>;

        return std::visit(
            Overload{
                [](std::conditional_t<C, std::add_const_t<T1>, T1>& p) -> R {
                  return p.weight;
                },
                [&](std::conditional_t<C, std::add_const_t<T2>, T2>& p) -> R {
                  return result.weightsOfStates.at(p.index());
                }},
            std::get<0>(cmp));
      };

      ////////////////////////
      // The Core Algorithm
      ////////////////////////

      if (haveMaterial) {
        if (haveMeasurement) {
          applyMultipleScattering(state, stepper, preUpdate);
        } else {
          applyMultipleScattering(state, stepper, fullUpdate);
        }
      }

      // If we have material, compute the bound state, split components
      // and reduce them to the desired size
      if (haveMaterial) {
        detail::extractComponents(
            state, stepper, result.parentTips,
            detail::ComponentSplitter{m_cfg.bethe_heitler_approx,
                                      m_cfg.weightCutoff},
            m_cfg.doCovTransport, componentCache);

        reduceComponents(stepper, surface, componentCache);
      }
      // If we have no material, just compute the bound state
      else {
        detail::extractComponents(state, stepper, result.parentTips,
                                  detail::ComponentForwarder{},
                                  m_cfg.doCovTransport, componentCache);
      }

      // If we have a measuerement, perform a kalman update
      if (haveMeasurement) {
        result.result = kalmanUpdate(state, found_source_link->second, result,
                                     componentCache);

        detail::computePosteriorWeights(componentCache, mapToProxyAndWeight);
      }

      // In every case, remove low weight components
      removeLowWeightComponents(componentCache, mapToWeight);

      // In the end, update the stepper
      auto updateRes = haveMeasurement
                           ? updateStepper(state, stepper, componentCache,
                                           mapProxyToWeightParsCov)
                           : updateStepper(state, stepper, componentCache,
                                           mapParsCacheToWeightParsCov);

      if (!updateRes.ok()) {
        SET_ERROR_AND_RETURN_OR_ABORT_ACTOR(updateRes.error());
      }
      result.parentTips = *updateRes;

      if (haveMeasurement) {
        result.currentTips = result.parentTips;
      }

      // If we only done preUpdate before, now do postUpdate
      if (haveMaterial && haveMeasurement) {
        applyMultipleScattering(state, stepper, postUpdate);
      }
    }

    throw_assert(detail::weightsAreNormalized(
                     stepper.constComponentIterable(state.stepping),
                     [](const auto& cmp) { return cmp.weight(); }),
                 "at the end not normalized");
  }

  template <typename stepper_t>
  void reduceComponents(const stepper_t& stepper, const Surface& surface,
                        std::vector<ComponentCache>& cmps) const {
    // Final component number
    const auto final_cmp_number = std::min(
        static_cast<std::size_t>(stepper.maxComponents), m_cfg.maxComponents);

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

        detail::reduceWithKLDistance(cmps, final_cmp_number,
                                     ParametersCacheProjector{}, angle_desc);
      } break;
      default: {
        detail::reduceWithKLDistance(cmps, final_cmp_number,
                                     ParametersCacheProjector{});
      }
    }
  }

  template <typename propagator_state_t, typename stepper_t>
  void removeMissedComponents(propagator_state_t& state,
                              const stepper_t& stepper,
                              std::vector<std::size_t>& current_tips) const {
    throw_assert(
        stepper.numberComponents(state.stepping) == current_tips.size(),
        "size mismatch");
    auto components = stepper.componentIterable(state.stepping);

    // 1) Compute the summed momentum and weight of the lost components
    double sumW_loss = 0.0;
    double sumWeightedQOverP_loss = 0.0;
    double initialQOverP = 0.0;

    for (const auto cmp : components) {
      if (cmp.status() != Intersection3D::Status::onSurface) {
        sumW_loss += cmp.weight();
        sumWeightedQOverP_loss += cmp.weight() * cmp.pars()[eFreeQOverP];
      }

      initialQOverP += cmp.weight() * cmp.pars()[eFreeQOverP];
    }

    throw_assert(
        sumW_loss < 1.0, "sumW_loss is 1, components:\n"
                             << [&]() {
                                  std::stringstream ss;
                                  for (const auto cmp : components) {
                                    ss << "cmp: w=" << cmp.weight()
                                       << ", s=" << cmp.status() << "\n";
                                  }
                                  return ss.str();
                                }());

    // 2) Adjust the momentum of the remaining components AND update the
    // current_tips vector
    double checkWeightSum = 0.0;
    double checkQOverPSum = 0.0;
    std::vector<std::size_t> new_tips;

    auto cmp_it = components.begin();
    auto tip_it = current_tips.begin();

    for (; tip_it != current_tips.end(); ++cmp_it, ++tip_it) {
      if ((*cmp_it).status() == Intersection3D::Status::onSurface) {
        auto& weight = (*cmp_it).weight();
        auto& qop = (*cmp_it).pars()[eFreeQOverP];

        weight /= (1.0 - sumW_loss);
        qop = qop * (1.0 - sumW_loss) + sumWeightedQOverP_loss;

        checkWeightSum += weight;
        checkQOverPSum += weight * qop;

        new_tips.push_back(*tip_it);
      }
    }

    current_tips = new_tips;

    // 3) Remove components
    stepper.removeMissedComponents(state.stepping);

    // 4) Some checks
    throw_assert(
        std::abs(checkQOverPSum - initialQOverP) < 1.e-4,
        "momentum mismatch, initial: " << std::setprecision(8) << initialQOverP
                                       << ", final: " << checkQOverPSum);

    throw_assert(
        std::abs(checkWeightSum - 1.0) < s_normalizationTolerance,
        "must sum up to 1 but is " << std::setprecision(8) << checkWeightSum);

    throw_assert(detail::weightsAreNormalized(
                     stepper.constComponentIterable(state.stepping),
                     [](const auto& cmp) { return cmp.weight(); }),
                 "not normalized");

    throw_assert(
        stepper.numberComponents(state.stepping) == current_tips.size(),
        "size mismatch");
  }

  /// Remove components with low weights and renormalize.
  /// TODO This function does not expect normalized components, but this
  /// could be redundant work...
  template <typename weight_projector_t>
  void removeLowWeightComponents(std::vector<ComponentCache>& cmps,
                                 const weight_projector_t& proj) const {
    detail::normalizeWeights(cmps, proj);

    cmps.erase(std::remove_if(
                   cmps.begin(), cmps.end(),
                   [&](auto& cmp) { return proj(cmp) < m_cfg.weightCutoff; }),
               cmps.end());

    detail::normalizeWeights(cmps, proj);
  }

  /// Function that updates the stepper with the component Cache
  /// @note Components with weight less than the weight-cutoff are ignored and not
  /// added to the stepper. The lost momentum is not compensated at the
  /// moment, but the components are reweighted
  template <typename propagator_state_t, typename stepper_t,
            typename projector_t>
  Result<std::vector<size_t>> updateStepper(
      propagator_state_t& state, const stepper_t& stepper,
      const std::vector<ComponentCache>& componentCache,
      const projector_t& proj) const {
    const auto& surface = *state.navigation.currentSurface;

    // We collect new tips in the loop
    std::vector<size_t> new_parent_tips;

    // Clear components before adding new ones
    stepper.clearComponents(state.stepping);

    // Finally loop over components
    for (const auto& [variant, meta] : componentCache) {
      const auto& [weight, pars, cov] = proj(variant);

      // Keep track of the indices
      std::visit(Overload{[&, &meta = meta](
                              const detail::GsfComponentParameterCache&) {
                            new_parent_tips.push_back(meta.parentIndex);
                          },
                          [&](const MultiTrajectory::TrackStateProxy& proxy) {
                            new_parent_tips.push_back(proxy.index());
                          }},
                 variant);

      // Add the component to the stepper
      const BoundTrackParameters bound(surface.getSharedPtr(), pars, cov);

      auto res = stepper.addComponent(state.stepping, std::move(bound), weight);

      if (!res.ok()) {
        return res.error();
      }

      auto& cmp = *res;
      cmp.jacobian() = meta.jacobian;
      cmp.jacToGlobal() = meta.jacToGlobal;
      cmp.pathAccumulated() = meta.pathLength;
      cmp.derivative() = meta.derivative;
      cmp.jacTransport() = meta.jacTransport;
    }

    return new_parent_tips;
  }

  template <typename propagator_state_t>
  Result<void> kalmanUpdate(const propagator_state_t& state,
                            const SourceLink& source_link, result_type& result,
                            std::vector<ComponentCache>& components) const {
    const auto& logger = state.options.logger;
    const auto& surface = *state.navigation.currentSurface;
    result.currentTips.clear();

    // Boolean flag, so that not every component increases the
    // result.measurementStates counter
    bool counted_as_measurement_state = false;

    for (auto& [variant, meta] : components) {
      // Create new track state
      result.currentTips.push_back(result.fittedStates.addTrackState(
          TrackStatePropMask::All, meta.parentIndex));

      const auto [weight, pars, cov] =
          std::get<detail::GsfComponentParameterCache>(variant);

      variant = result.fittedStates.getTrackState(result.currentTips.back());
      auto& trackProxy = std::get<TrackProxy>(variant);

      // Set track parameters
      trackProxy.predicted() = std::move(pars);

      if (cov) {
        trackProxy.predictedCovariance() = std::move(*cov);
      }
      result.weightsOfStates[result.currentTips.back()] = weight;

      // Set surface
      trackProxy.setReferenceSurface(surface.getSharedPtr());

      // assign the source link to the track state
      trackProxy.setUncalibrated(source_link);

      trackProxy.jacobian() = std::move(meta.jacobian);
      trackProxy.pathLength() = std::move(meta.pathLength);

      // We have predicted parameters, so calibrate the uncalibrated
      m_cfg.extensions.calibrator(state.geoContext, trackProxy);

      // Get and set the type flags
      trackProxy.typeFlags().set(TrackStateFlag::ParameterFlag);
      if (surface.surfaceMaterial() != nullptr) {
        trackProxy.typeFlags().set(TrackStateFlag::MaterialFlag);
      }

      // Do Kalman update
      if (not m_cfg.extensions.outlierFinder(trackProxy)) {
        // Perform update
        auto updateRes = m_cfg.extensions.updater(state.geoContext, trackProxy,
                                   state.stepping.navDir, logger);

        if (!updateRes.ok()) {
          ACTS_ERROR("Update step failed: " << updateRes.error());
          RETURN_ERROR_OR_ABORT_ACTOR(updateRes.error());
        }

        trackProxy.typeFlags().set(TrackStateFlag::MeasurementFlag);

        // We count the state with measurement TODO does this metric make
        // sense for a GSF?
        if (!counted_as_measurement_state) {
          ++result.measurementStates;
          counted_as_measurement_state = true;
        }
      } else {
        // TODO ACTS_VERBOSE message
        throw std::runtime_error("outlier handling not yet implemented fully");

        // Set the outlier type flag
        trackProxy.typeFlags().set(TrackStateFlag::OutlierFlag);
        trackProxy.data().ifiltered = trackProxy.data().ipredicted;
      }
    }

    ++result.processedStates;

    return Acts::Result<void>::success();
  }

  /// Apply the multipe scattering to the state
  template <typename propagator_state_t, typename stepper_t>
  void applyMultipleScattering(
      propagator_state_t& state, const stepper_t& stepper,
      const MaterialUpdateStage& updateStage = fullUpdate) const {
    const auto& logger = state.options.logger;
    const auto& surface = *state.navigation.currentSurface;

    for (auto cmp : stepper.componentIterable(state.stepping)) {
      auto singleState = cmp.singleState(state);
      const auto &singleStepper = cmp.singleStepper(stepper);

      detail::PointwiseMaterialInteraction interaction(&surface, singleState,
                                                       singleStepper);

      if (interaction.evaluateMaterialSlab(singleState, updateStage)) {
        constexpr bool doMultipleScattering = true;
        constexpr bool doEnergyLoss = false;
        interaction.evaluatePointwiseMaterialInteraction(doMultipleScattering,
                                                         doEnergyLoss);

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
}  // namespace Acts

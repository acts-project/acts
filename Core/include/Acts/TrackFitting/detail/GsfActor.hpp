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
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/BetheHeitlerApprox.hpp"
#include "Acts/TrackFitting/detail/GsfSmoothing.hpp"
#include "Acts/TrackFitting/detail/GsfUtils.hpp"
#include "Acts/TrackFitting/detail/KLMixtureReduction.hpp"
#include "Acts/TrackFitting/detail/kalman_update_helpers.hpp"
#include "Acts/Utilities/Overload.hpp"

#include <ios>
#include <map>
#include <numeric>

namespace Acts {
namespace detail {

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

    /// Wether to transport covariance
    bool doCovTransport = true;

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
    KalmanFitterExtensions extensions;
  } m_cfg;

  /// Broadcast Cache Type
  using TrackProxy = typename MultiTrajectory::TrackStateProxy;
  using ComponentCache = std::tuple<detail::GsfComponentParameterCache,
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
    // Set error or abort utility
    auto set_error_or_abort = [&](auto error) {
      if (m_cfg.abortOnError) {
        std::abort();
      } else {
        result.result = error;
      }
    };

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

    // Workaround to initialize e.g. MultiTrajectory in backward mode
    if (!result.haveInitializedResult && m_cfg.resultInitializer) {
      m_cfg.resultInitializer(result, logger);
      result.haveInitializedResult = true;
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
          !m_cfg.disableAllMaterialHandling;
      const bool haveMeasurement =
          found_source_link != m_cfg.inputMeasurements.end();

      ACTS_VERBOSE(std::boolalpha << "haveMaterial " << haveMaterial
                                  << ", haveMeasurement: " << haveMeasurement);

      // Collect holes in result, but no holes before first measuerement
      if (result.measurementStates > 0 && !haveMeasurement &&
          surface.associatedDetectorElement()) {
        ACTS_VERBOSE("Detected hole on " << surface.geometryId());
        result.missedActiveSurfaces.push_back(&surface);
        ++result.measurementHoles;
      }

      // Early return if nothing happens
      if (not haveMaterial && not haveMeasurement) {
        return;
      }

      // Projectors
      //       auto mapToProxyAndWeight = [&](auto& cmp) {
      //         auto& proxy = std::get<TrackProxy>(std::get<0>(cmp));
      //         return std::tie(proxy,
      //         result.weightsOfStates.at(proxy.index()));
      //       };
      //
      //       auto mapParsCacheToWeightParsCov = [&](auto& variant) ->
      //       decltype(auto) {
      //         return std::get<detail::GsfComponentParameterCache>(variant);
      //       };
      //
      //       auto mapProxyToWeightParsCov = [&](auto& variant) {
      //         auto& proxy = std::get<TrackProxy>(variant);
      //         return
      //         std::make_tuple(result.weightsOfStates.at(proxy.index()),
      //                                proxy.filtered(),
      //                                proxy.filteredCovariance());
      //       };
      //
      //       auto mapToWeight = [&](auto& cmp) -> decltype(auto) {
      //         constexpr bool C =
      //             std::is_const_v<std::remove_reference_t<decltype(cmp)>>;
      //
      //         using T1 = detail::GsfComponentParameterCache;
      //         using T2 = MultiTrajectory::TrackStateProxy;
      //         using R = std::conditional_t<C, double, double&>;
      //
      //         return std::visit(
      //             Overload{
      //                 [](std::conditional_t<C, std::add_const_t<T1>, T1>& p)
      //                 -> R {
      //                   return p.weight;
      //                 },
      //                 [&](std::conditional_t<C, std::add_const_t<T2>, T2>& p)
      //                 -> R {
      //                   return result.weightsOfStates.at(p.index());
      //                 }},
      //             std::get<0>(cmp));
      //       };

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

      // We do not need the component cache here, we can just update our stepper
      // state with the filtered components.
      // NOTE because of early return before we know that we have a measurement
      if (not haveMaterial) {
        auto update_stepper = [&](auto cmp, const auto& proxy, double) {
          cmp.pars() = MultiTrajectoryHelpers::freeFiltered(
              state.options.geoContext, proxy);
          cmp.cov() = proxy.filteredCovariance();
        };

        kalmanUpdate(state, stepper, result, found_source_link->second,
                     update_stepper);
      }

      if (haveMaterial) {
        // Create Cache
        std::vector<ComponentCache> componentCache;

        // Here we create the new with the filtered kalman result
        if (haveMeasurement) {
          auto update_cache = [&](auto cmp, const auto& proxy, double weight) {
            detail::GsfComponentMetaCache mcache;
            mcache.parentIndex = proxy.index();
            mcache.jacobian = cmp.jacobian();
            mcache.jacToGlobal = cmp.jacToGlobal();
            mcache.jacTransport = cmp.jacTransport();
            mcache.derivative = cmp.derivative();
            mcache.pathLength = cmp.pathAccumulated();

            BoundTrackParameters bound(surface.getSharedPtr(), proxy.filtered(),
                                       proxy.filteredCovariance());
            detail::create_new_components(state, bound, weight, mcache,
                                          m_cfg.bethe_heitler_approx,
                                          componentCache, m_cfg.weightCutoff);
          };

          kalmanUpdate(state, stepper, result, found_source_link->second,
                       update_cache);
        } 
        // Here we simply compute the bound state and then create the new components
        else {
          auto cmp_iterable = stepper.componentIterable(state.stepping);
          auto cmp_it = cmp_iterable.begin();
          auto idx_it = result.parentTips.begin();

          for (; idx_it != result.parentTips.end(); ++cmp_it, ++idx_it) {
            auto cmp = *cmp_it;

            detail::GsfComponentMetaCache mcache;
            mcache.parentIndex = *idx_it;
            mcache.jacobian = cmp.jacobian();
            mcache.jacToGlobal = cmp.jacToGlobal();
            mcache.jacTransport = cmp.jacTransport();
            mcache.derivative = cmp.derivative();
            mcache.pathLength = cmp.pathAccumulated();

            if (cmp.status() != Intersection3D::Status::onSurface) {
              ACTS_VERBOSE("Skip component which is not on surface");
              continue;
            }

            auto boundState = cmp.boundState(surface, m_cfg.doCovTransport);

            if (!boundState.ok()) {
              ACTS_ERROR(
                  "Failed to compute boundState: " << boundState.error());
              continue;
            }

            const auto& [bound, jac, pathLength] = boundState.value();

            detail::create_new_components(state, bound, cmp.weight(), mcache,
                                          m_cfg.bethe_heitler_approx,
                                          componentCache, m_cfg.weightCutoff);
          }
        }
        
        reduceComponents(stepper, surface, componentCache);
      }

      // If we have material, compute the bound state, split components
      // and reduce them to the desired size
      if (haveMaterial) {
        detail::extractComponents(
            state, stepper, result.parentTips,
            detail::ComponentSplitter{m_cfg.bethe_heitler_approx,
                                      m_cfg.weightCutoff},
            m_cfg.doCovTransport, componentCache);

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
        return set_error_or_abort(updateRes.error());
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

  template <typename propagator_state_t, typename stepper_t,
            typename component_handler_t>
  Result<void> kalmanUpdate(propagator_state_t& state, const stepper_t& stepper,
                            result_type& result, const SourceLink& source_link,
                            component_handler_t& handle_component) const {
    const auto& logger = state.options.logger;
    const auto& surface = *state.navigation.currentSurface;

    // Clear things that are overwritten soon
    result.currentTips.clear();

    // Boolean flag, to distinguish measurement and outlier states. This flag
    // is only modified by the valid-measurement-branch, so only if there
    // isn't any valid measuerement state, the flag stays false and the state
    // is thus counted as an outlier
    bool is_valid_measurement = false;

    // Do a for loop over the components and the parent indices
    auto cmp_iterable = stepper.componentIterable(state.stepping);
    auto cmp_it = cmp_iterable.begin();
    auto idx_it = result.parentTips.begin();

    for (; idx_it != result.parentTips.end(); ++cmp_it, ++idx_it) {
      auto cmp = *cmp_it;

      if (cmp.status() != Intersection3D::Status::onSurface) {
        ACTS_VERBOSE("Skip component which is not on surface");
        continue;
      }

      auto singleState = cmp.singleState(state);
      const auto& singleStepper = cmp.singleStepper(stepper);

      auto trackStateProxyRes = detail::handleMeasurement(
          singleState, singleStepper, m_cfg.extensions, surface, source_link,
          result.fittedStates, *idx_it);

      if (!trackStateProxyRes.ok()) {
        return trackStateProxyRes.error();
      }

      const auto& trackStateProxy = *trackStateProxyRes;
      result.currentTips.push_back(trackStateProxy.index());

      // If at least one component is no outlier, we consider the whole thing
      // as a measuerementState
      if (trackStateProxy.typeFlags().test(
              Acts::TrackStateFlag::MeasurementFlag)) {
        is_valid_measurement = true;
      }

      // Do whatever needs to be done next
      handle_component(cmp, trackStateProxy,
                       result.weightsOfStates.at(*idx_it));
    }

    // Do the statistics
    ++result.processedStates;

    if (is_valid_measurement) {
      ++result.measurementStates;
    }

    // Return sucess
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

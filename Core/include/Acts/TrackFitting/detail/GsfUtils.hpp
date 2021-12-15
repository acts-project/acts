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
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>
#include <numeric>

namespace Acts {

/// The tolerated difference to 1 to accept weights as normalized
/// TODO seems sometimes to fail for 1.e-8
constexpr static double s_normalizationTolerance = 1.e-4;

namespace detail {

template <typename component_range_t, typename projector_t,
          typename print_flag_t = std::false_type>
bool weightsAreNormalized(const component_range_t &cmps,
                          const projector_t &proj,
                          double tol = s_normalizationTolerance,
                          print_flag_t print_flag = print_flag_t{}) {
  double sum_of_weights = 0.0;

  for (const auto &cmp : cmps) {
    sum_of_weights += proj(cmp);
  }

  if (std::abs(sum_of_weights - 1.0) < tol) {
    return true;
  } else {
    if constexpr (print_flag) {
      std::cout << std::setprecision(10)
                << "diff from 1: " << std::abs(sum_of_weights - 1.0) << "\n";
    }

    return false;
  }
}

template <typename component_range_t, typename projector_t>
void normalizeWeights(component_range_t &cmps, const projector_t &proj) {
  double sum_of_weights = 0.0;

  for(auto it=cmps.begin(); it != cmps.end(); ++it) {
    sum_of_weights += proj(*it);
  }

  for(auto it=cmps.begin(); it != cmps.end(); ++it) {
    proj(*it) /= sum_of_weights;
  }
}

/// Stores meta information about the components
struct GsfComponentMetaCache {
  /// Where to find the parent component in the MultiTrajectory
  std::size_t parentIndex;

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
struct GsfComponentParameterCache {
  ActsScalar weight;
  BoundVector boundPars;
  std::optional<BoundSymMatrix> boundCov;
};

/// @brief A multi component state
using MultiComponentState =
    std::pair<std::shared_ptr<const Surface>,
              std::vector<std::tuple<ActsScalar, BoundVector,
                                     std::optional<BoundSymMatrix>>>>;

/// Functor which splits a component into multiple new components by creating a
/// gaussian mixture, and adds them to a cache with a customizable cache_maker
template <typename bethe_heitler_t>
struct ComponentSplitter {
  const bethe_heitler_t &betheHeitler;
  const double weightCutoff;

  ComponentSplitter(const bethe_heitler_t &bh, double cutoff)
      : betheHeitler(bh), weightCutoff(cutoff) {}

  template <typename propagator_state_t, typename component_cache_t>
  void operator()(const propagator_state_t &state,
                  const BoundTrackParameters &old_bound,
                  const double old_weight,
                  const GsfComponentMetaCache &metaCache,
                  std::vector<component_cache_t> &componentCaches) const {
    const auto &logger = state.options.logger;
    const auto &surface = *state.navigation.currentSurface;
    const auto p_prev = old_bound.absoluteMomentum();

    // Evaluate material slab
    auto slab =
        surface.surfaceMaterial()->materialSlab(
            old_bound.position(state.stepping.geoContext),
            state.stepping.navDir, MaterialUpdateStage::fullUpdate);

//     auto pathCorrection = surface.pathCorrection(
//         state.stepping.geoContext, old_bound.position(state.stepping.geoContext),
//         old_bound.unitDirection());
//     slab.scaleThickness(pathCorrection);

    // Get the mixture
    const auto mixture = betheHeitler->mixture(slab.thicknessInX0());

    // Create all possible new components
    for (const auto &gaussian : mixture) {
      // Here we combine the new child weight with the parent weight.
      // However, this must be later re-adjusted
      const auto new_weight = gaussian.weight * old_weight;

      if (new_weight < weightCutoff) {
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
        if (state.stepping.navDir == NavigationDirection::forward)
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
          if (state.stepping.navDir == NavigationDirection::forward) {
            const auto f = 1. / (p_prev * gaussian.mean);
            return f * f * gaussian.var;
          } else {
            return gaussian.var / (p_prev * p_prev);
          }
        }();

        (*new_cov)(eBoundQOverP, eBoundQOverP) += varInvP;
        throw_assert(std::isfinite((*new_cov)(eBoundQOverP, eBoundQOverP)), "cov not finite, varInvP=" << varInvP << ", p_prev=" << p_prev << ", gaussian.mean=" << gaussian.mean << ", gaussian.var=" << gaussian.var);
      }

      // Set the remaining things and push to vector
      componentCaches.push_back(
          {GsfComponentParameterCache{new_weight, new_pars, new_cov},
           metaCache});
    }
  }
};

/// Functor whicht forwards components to the component cache and does not
/// splitting therefore
struct ComponentForwarder {
  template <typename propagator_state_t, typename component_cache_t,
            typename meta_cache_t>
  void operator()(const propagator_state_t &,
                  const BoundTrackParameters &old_bound,
                  const double old_weight, const meta_cache_t &metaCache,
                  std::vector<component_cache_t> &componentCaches) const {
    componentCaches.push_back(
        {GsfComponentParameterCache{old_weight, old_bound.parameters(),
                                    old_bound.covariance()},
         metaCache});
  }
};

/// @brief Expands all existing components to new components by using a
/// gaussian-mixture approximation for the Bethe-Heitler distribution.
/// @return a std::vector with all new components (parent tip, weight,
/// parameters, covariance)
template <typename propagator_state_t, typename stepper_t, typename component_t,
          typename component_processor_t>
void extractComponents(propagator_state_t &state, const stepper_t &stepper,
                       const std::vector<std::size_t> &parentTrajectoryIdxs,
                       const component_processor_t &componentProcessor,
                       const bool doCovTransport,
                       std::vector<component_t> &componentCache) {
  // Some shortcuts
  auto &stepping = state.stepping;
  const auto &logger = state.options.logger;
  const auto &surface = *state.navigation.currentSurface;

  // Approximate bethe-heitler distribution as gaussian mixture
  std::size_t i = 0;
  for (auto old_cmp : stepper.componentIterable(stepping)) {
    if (old_cmp.status() != Intersection3D::Status::onSurface) {
      ACTS_VERBOSE("Skip component which is not on surface");
      continue;
    }

    auto boundState = old_cmp.boundState(surface, doCovTransport);

    if (!boundState.ok()) {
      ACTS_ERROR("Failed to compute boundState: " << boundState.error());
      continue;
    }

    const auto &[old_bound, jac, pathLength] = boundState.value();

    detail::GsfComponentMetaCache metaCache{
        parentTrajectoryIdxs[i++], jac,
        old_cmp.jacToGlobal(),     old_cmp.jacTransport(),
        old_cmp.derivative(),      pathLength};

    componentProcessor(state, old_bound, old_cmp.weight(), metaCache,
                       componentCache);
  }
}

/// @brief Function object which maps a value to itself by perfect forwarding
/// TODO Replace Identity with std::identity on C++20
struct Identity {
  template <typename T>
  auto operator()(T &&v) const {
    return std::forward<T>(v);
  }
};

/// @brief Reweight the components according to `R. Frühwirth, "Track fitting
/// with non-Gaussian noise"`. See also the implementation in Athena at
/// PosteriorWeightsCalculator.cxx
/// Expects that the projector maps the component to something like a
/// std::pair< trackProxy&, double& > so that it can be extracted with std::get
/// @note The weights are not renormalized!
template <typename component_t, typename projector_t>
void computePosteriorWeights(std::vector<component_t> &cmps,
                             const projector_t &proj) {
  // Helper Function to compute detR
  auto computeDetR = [](const auto &trackState) -> ActsScalar {
    const auto predictedCovariance = trackState.predictedCovariance();

    return visit_measurement(
        trackState.calibrated(), trackState.calibratedCovariance(),
        trackState.calibratedSize(),
        [&](const auto calibrated, const auto calibratedCovariance) {
          constexpr size_t kMeasurementSize =
              decltype(calibrated)::RowsAtCompileTime;
          const auto H =
              trackState.projector()
                  .template topLeftCorner<kMeasurementSize, eBoundSize>()
                  .eval();

          return (H * predictedCovariance * H.transpose() +
                  calibratedCovariance)
              .determinant();
        });
  };

  // Find minChi2, this can be used to factor some things later in the
  // exponentiation
  const auto minChi2 =
      std::get<0>(proj(*std::min_element(cmps.begin(), cmps.end(),
                                         [&](const auto &a, const auto &b) {
                                           return std::get<0>(proj(a)).chi2() <
                                                  std::get<0>(proj(b)).chi2();
                                         })))
          .chi2();

  // Compute new weights and reweight
  double sumOfWeights = 0.;

  for (auto &cmp : cmps) {
    const double chi2 = std::get<0>(proj(cmp)).chi2() - minChi2;
    const double detR = computeDetR(std::get<0>(proj(cmp)));

    if (!std::isfinite(chi2) || !std::isfinite(detR)) {
      sumOfWeights += std::get<1>(proj(cmp));
    } else {
      const auto factor = std::sqrt(1. / detR) * std::exp(-0.5 * chi2);
      std::get<1>(proj(cmp)) *= factor;
    }
  }
}

/// Enumeration type used in extractMultiComponentStates(...)
enum class StatesType { ePredicted, eFiltered, eSmoothed };

inline std::ostream &operator<<(std::ostream &os, StatesType type) {
  constexpr static std::array names = {"predicted", "filtered", "smoothed"};
  os << names[static_cast<int>(type)];
  return os;
}

/// @brief Extracts a MultiComponentState from a MultiTrajectory and a given list of indices
auto extractMultiComponentState(const MultiTrajectory &traj,
                                const std::vector<size_t> &tips,
                                const std::map<size_t, ActsScalar> &weights,
                                StatesType type)
    -> MultiComponentBoundTrackParameters<SinglyCharged>;

}  // namespace detail

}  // namespace Acts

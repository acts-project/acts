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

  for (auto it = cmps.begin(); it != cmps.end(); ++it) {
    sum_of_weights += proj(*it);
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

  // we need decltype(auto) here to support proxy-types with reference
  // semantics, otherwise there is a `cannot bind ... to ...` error
  for (auto it = cmps.begin(); it != cmps.end(); ++it) {
    decltype(auto) cmp = *it;
    sum_of_weights += proj(cmp);
  }

  for (auto it = cmps.begin(); it != cmps.end(); ++it) {
    decltype(auto) cmp = *it;
    proj(cmp) /= sum_of_weights;
  }
}

// A class that prints information about the state on construction and
// destruction, it also contains some assertions in the constructor and
// destructor. It can be removed without change of behaviour, since it only
// holds const references
template <typename propagator_state_t, typename stepper_t>
class ScopedGsfInfoPrinterAndChecker {
  const propagator_state_t &m_state;
  const stepper_t &m_stepper;
  double m_p_initial;

  const auto &logger() const { return m_state.options.logger(); }

  void print_component_stats() const {
    std::size_t i = 0;
    for (auto cmp : m_stepper.constComponentIterable(m_state.stepping)) {
      auto getVector = [&](auto idx) {
        return cmp.pars().template segment<3>(idx).transpose();
      };
      ACTS_VERBOSE("  #" << i++ << " pos: " << getVector(eFreePos0) << ", dir: "
                         << getVector(eFreeDir0) << ", weight: " << cmp.weight()
                         << ", status: " << cmp.status()
                         << ", qop: " << cmp.pars()[eFreeQOverP]);
    }
  }

  void checks(const std::string_view &where) const {
    const auto cmps = m_stepper.constComponentIterable(m_state.stepping);
    throw_assert(detail::weightsAreNormalized(
                     cmps, [](const auto &cmp) { return cmp.weight(); }),
                 "not normalized at " << where);

    throw_assert(
        std::all_of(cmps.begin(), cmps.end(),
                    [](auto cmp) { return std::isfinite(cmp.weight()); }),
        "some weights are not finite at " << where);
  }

 public:
  ScopedGsfInfoPrinterAndChecker(const propagator_state_t &state,
                                const stepper_t &stepper)
      : m_state(state),
        m_stepper(stepper),
        m_p_initial(stepper.momentum(state.stepping)) {
    // Some initial printing
    checks("start");
    ACTS_VERBOSE("Gsf step "
                 << state.stepping.steps << " at mean position "
                 << stepper.position(state.stepping).transpose()
                 << " with direction "
                 << stepper.direction(state.stepping).transpose()
                 << " and momentum " << stepper.momentum(state.stepping)
                 << " and charge " << stepper.charge(state.stepping));
    ACTS_VERBOSE("Propagation is in "
                 << (state.stepping.navDir == forward ? "forward" : "backward")
                 << " mode");
    print_component_stats();
  }

  ~ScopedGsfInfoPrinterAndChecker() {
    if (m_state.navigation.currentSurface) {
      const auto p_final = m_stepper.momentum(m_state.stepping);
      ACTS_VERBOSE("Component status at end of step:");
      print_component_stats();
      ACTS_VERBOSE("Delta Momentum = " << std::setprecision(5)
                                       << p_final - m_p_initial);
    }
    checks("end");
  }
};

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

/// Reweight the components according to `R. Frühwirth, "Track fitting
/// with non-Gaussian noise"`. See also the implementation in Athena at
/// PosteriorWeightsCalculator.cxx
/// @note The weights are not renormalized!
void computePosteriorWeights(const MultiTrajectory &mt,
                             const std::vector<std::size_t> &tips,
                             std::map<std::size_t, double> &weights);

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

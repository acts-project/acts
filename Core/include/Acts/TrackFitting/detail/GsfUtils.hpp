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
  std::size_t m_missedCount;

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

    // If all components are missed, their weights have been reset to zero.
    // In this case the weights might not be normalized and not even be
    // finite due to a division by zero.
    if (m_stepper.numberComponents(m_state.stepping) > m_missedCount) {
      throw_assert(detail::weightsAreNormalized(
                       cmps, [](const auto &cmp) { return cmp.weight(); }),
                   "not normalized at " << where);

      throw_assert(
          std::all_of(cmps.begin(), cmps.end(),
                      [](auto cmp) { return std::isfinite(cmp.weight()); }),
          "some weights are not finite at " << where);
    }
  }

 public:
  ScopedGsfInfoPrinterAndChecker(const propagator_state_t &state,
                                 const stepper_t &stepper,
                                 std::size_t missedCount)
      : m_state(state),
        m_stepper(stepper),
        m_p_initial(stepper.momentum(state.stepping)),
        m_missedCount(missedCount) {
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
                 << (state.stepping.navDir == NavigationDirection::Forward
                         ? "forward"
                         : "backward")
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

ActsScalar calculateDeterminant(
    const double *fullCalibrated, const double *fullCalibratedCovariance,
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     true>::Covariance predictedCovariance,
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax, true>::Projector
        projector,
    unsigned int calibratedSize);

/// Reweight the components according to `R. Frühwirth, "Track fitting
/// with non-Gaussian noise"`. See also the implementation in Athena at
/// PosteriorWeightsCalculator.cxx
/// @note The weights are not renormalized!
template <typename D>
void computePosteriorWeights(
    const MultiTrajectory<D> &mt,
    const std::vector<MultiTrajectoryTraits::IndexType> &tips,
    std::map<MultiTrajectoryTraits::IndexType, double> &weights) {
  // Helper Function to compute detR

  // Find minChi2, this can be used to factor some things later in the
  // exponentiation
  const auto minChi2 =
      mt.getTrackState(*std::min_element(tips.begin(), tips.end(),
                                         [&](const auto &a, const auto &b) {
                                           return mt.getTrackState(a).chi2() <
                                                  mt.getTrackState(b).chi2();
                                         }))
          .chi2();

  // Loop over the tips and compute new weights
  for (auto tip : tips) {
    const auto state = mt.getTrackState(tip);
    const double chi2 = state.chi2() - minChi2;
    const double detR = calculateDeterminant(
        // This abuses an incorrectly sized vector / matrix to access the
        // data pointer! This works (don't use the matrix as is!), but be
        // careful!
        state.template calibrated<MultiTrajectoryTraits::MeasurementSizeMax>()
            .data(),
        state
            .template calibratedCovariance<
                MultiTrajectoryTraits::MeasurementSizeMax>()
            .data(),
        state.predictedCovariance(), state.projector(), state.calibratedSize());

    const auto factor = std::sqrt(1. / detR) * std::exp(-0.5 * chi2);

    // If something is not finite here, just leave the weight as it is
    if (std::isfinite(factor)) {
      weights.at(tip) *= factor;
    }
  }
}

}  // namespace detail
}  // namespace Acts

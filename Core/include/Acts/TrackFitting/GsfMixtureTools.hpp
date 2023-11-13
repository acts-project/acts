// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/detail/MergeGaussianMixture.hpp"

namespace Acts {

/// A component in a Gaussian mixture
struct GsfComponent {
  ActsScalar weight = 0;
  BoundVector boundPars = BoundVector::Zero();
  BoundSquareMatrix boundCov = BoundSquareMatrix::Identity();
};

/// Very simple mixture reduction method: Just removes the components with the
/// smallest weight until the required number of components is reached
/// @param cmpCache the component collection
/// 
/// @param maxCmpsAfterMerge the number of components we want to reach
/// @param surface the surface type on which the components are (unused here)
void reduceMixtureLargestWeights(std::vector<Acts::GsfComponent> &cmpCache,
                                 std::size_t maxCmpsAfterMerge,
                                 const Surface &surface);

/// Greedy component reduction algorithm. Reduces the components with the
/// minimal symmetric KL-distance (applied only to the q/p-dimension) until the
/// required number of components is reached.
/// 
/// @param cmpCache the component collection
/// @param maxCmpsAfterMerge the number of components we want to reach
/// @param surface the surface type on which the components are
void reduceMixtureWithKLDistance(std::vector<Acts::GsfComponent> &cmpCache,
                                 std::size_t maxCmpsAfterMerge,
                                 const Surface &surface);

/// @enum ComponentMergeMethod
///
/// Available reduction methods for the reduction of a Gaussian mixture
enum class ComponentMergeMethod { eMean, eMaxWeight };

/// Reduce Gaussian mixture
///
/// @param mixture The mixture iterable
/// @param surface The surface, on which the bound state is
/// @param method How to reduce the mixture
/// @param projector How to project a element of the iterable to something
/// like a std::tuple< double, BoundVector, BoundMatrix >
///
/// @return parameters and covariance as std::tuple< BoundVector, BoundMatrix >
template <typename mixture_t, typename projector_t = Acts::Identity>
auto mergeGaussianMixture(const mixture_t &mixture, const Surface &surface,
                          ComponentMergeMethod method,
                          projector_t &&projector = projector_t{}) {
  using T = std::tuple<Acts::BoundVector, Acts::BoundSquareMatrix>;
  using R = Acts::Result<T>;

  return detail::angleDescriptionSwitch(surface, [&](const auto &desc) {
    BoundVector parameters = BoundVector::Zero();
    switch (method) {
      case ComponentMergeMethod::eMean: {
        parameters = detail::computeMixtureMean(mixture, projector, desc);
      } break;
      case ComponentMergeMethod::eMaxWeight: {
        const auto maxWeightIt = std::max_element(
            mixture.begin(), mixture.end(), [&](const auto &a, const auto &b) {
              return std::get<0>(projector(a)) < std::get<0>(projector(b));
            });
        parameters = std::get<1>(projector(*maxWeightIt));
      }
    }

    // MARK: fpeMaskBegin(FLTUND, 1, #2347)
    const auto covariance =
        computeMixtureCovariance(mixture, parameters, projector, desc);
    // MARK: fpeMaskEnd(FLTUND)

    return R{T{parameters, covariance}};
  });
}

}  // namespace Acts

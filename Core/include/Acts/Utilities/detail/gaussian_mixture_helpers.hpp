// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/Identity.hpp"

#include <cmath>
#include <tuple>

namespace Acts {
namespace detail {

/// Namespace containing Angle descriptions for the combineBoundGaussianMixture
/// function
namespace AngleDescription {

template <BoundIndices Idx>
struct CyclicAngle {
  constexpr static BoundIndices idx = Idx;
  constexpr static double constant = 1.0;
};

template <BoundIndices Idx>
struct CyclicRadiusAngle {
  constexpr static BoundIndices idx = Idx;
  double constant = 1.0;  // the radius
};

using Default = std::tuple<CyclicAngle<eBoundPhi>>;
using Cylinder =
    std::tuple<CyclicRadiusAngle<eBoundLoc0>, CyclicAngle<eBoundPhi>>;
}  // namespace AngleDescription

/// @brief Combine multiple components into one representative track state
/// object. The function takes iterators to allow for arbitrary ranges to be
/// combined
/// @tparam component_iterator_t An iterator of a range of components
/// @tparam projector_t A projector, which maps the component to a std::tuple<
/// ActsScalar, BoundVector, std::optional< BoundSymMatrix > >
/// @tparam angle_desc_t A angle description object which defines the cyclic
/// angles in the bound parameters
template <typename component_iterator_t, typename projector_t = Identity,
          typename angle_desc_t = AngleDescription::Default>
auto combineBoundGaussianMixture(
    const component_iterator_t begin, const component_iterator_t end,
    projector_t &&projector = projector_t{},
    const angle_desc_t &angle_desc = angle_desc_t{}) {
  throw_assert(std::distance(begin, end) > 0, "empty component range");

  using ret_type = std::tuple<BoundVector, std::optional<BoundSymMatrix>>;
  BoundVector mean = BoundVector::Zero();
  BoundSymMatrix cov1 = BoundSymMatrix::Zero();
  BoundSymMatrix cov2 = BoundSymMatrix::Zero();
  double sumOfWeights{0.0};

  const auto &[begin_weight, begin_pars, begin_cov] = projector(*begin);

  if (std::distance(begin, end) == 1) {
    return ret_type{begin_pars, *begin_cov};
  }

  // clang-format off
  // x = \sum_{l} w_l * x_l
  // C = \sum_{l} w_l * C_l + \sum_{l} \sum_{m>l} w_l * w_m * (x_l - x_m)(x_l - x_m)^T
  // clang-format on
  for (auto l = begin; l != end; ++l) {
    const auto &[weight_l, pars_l, cov_l] = projector(*l);
    throw_assert(cov_l, "we require a covariance here");

    sumOfWeights += weight_l;
    mean += weight_l * pars_l;
    cov1 += weight_l * *cov_l;

    // Avoid problems with cyclic coordinates. The indices for these are taken
    // from the angle_description_t template parameter, and applied with the
    // following lambda.
    auto handleCyclicCoor = [&ref = begin_pars, &pars = pars_l,
                             &weight = weight_l, &mean = mean](auto desc) {
      const auto delta = (ref[desc.idx] - pars[desc.idx]) / desc.constant;

      if (delta > M_PI) {
        mean[desc.idx] += (2 * M_PI) * weight * desc.constant;
      } else if (delta < -M_PI) {
        mean[desc.idx] -= (2 * M_PI) * weight * desc.constant;
      }
    };

    std::apply([&](auto... idx_desc) { (handleCyclicCoor(idx_desc), ...); },
               angle_desc);

    // For covariance we must loop over all other following components
    for (auto m = std::next(l); m != end; ++m) {
      const auto &[weight_m, pars_m, cov_m] = projector(*m);
      throw_assert(cov_m, "we require a covariance here");

      const BoundVector diff = pars_l - pars_m;
      cov2 += weight_l * weight_m * diff * diff.transpose();
    }
  }

  return ret_type{mean / sumOfWeights,
                  cov1 / sumOfWeights + cov2 / (sumOfWeights * sumOfWeights)};
}

}  // namespace detail
}  // namespace Acts

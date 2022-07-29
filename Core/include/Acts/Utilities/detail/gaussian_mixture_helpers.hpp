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
#include "Acts/Utilities/detail/periodic.hpp"

#include <cmath>
#include <optional>
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
using Disk = std::tuple<CyclicAngle<eBoundLoc1>, CyclicAngle<eBoundPhi>>;
}  // namespace AngleDescription

/// @brief Combine multiple components into one representative track state
/// object. The function takes iterators to allow for arbitrary ranges to be
/// combined. The dimension of the vectors is infeared from the inputs.
///
/// @note If one component does not contain a covariance, no covariance is
/// computed.
///
/// @note The correct mean and variances for cyclic coordnates or spherical
/// coordinates (theta, phi) must generally be computed using a special circular
/// mean or in cartesian coordinates. This implements a approximation, which
/// only works well for close components.
///
/// @tparam component_iterator_t An iterator of a range of components
/// @tparam projector_t A projector, which maps the component to a
/// std::tuple< weight, mean, std::optional< cov > >
/// @tparam angle_desc_t A angle description object which defines the cyclic
/// angles in the bound parameters
template <typename component_iterator_t, typename projector_t = Identity,
          typename angle_desc_t = AngleDescription::Default>
auto combineGaussianMixture(const component_iterator_t begin,
                            const component_iterator_t end,
                            projector_t &&projector = projector_t{},
                            const angle_desc_t &angleDesc = angle_desc_t{}) {
  // Extract the first component
  const auto &[beginWeight, beginPars, beginCov] = projector(*begin);

  // Assert type properties
  using ParsType = std::decay_t<decltype(beginPars)>;
  using CovType = std::decay_t<decltype(*beginCov)>;
  using WeightType = std::decay_t<decltype(beginWeight)>;

  constexpr int D = ParsType::RowsAtCompileTime;
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(ParsType);
  EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(CovType, D, D);
  static_assert(std::is_floating_point_v<WeightType>);

  std::apply(
      [&](auto... d) { static_assert((std::less<int>{}(d.idx, D) && ...)); },
      angleDesc);

  // Define the return type
  using RetType = std::tuple<ActsVector<D>, std::optional<ActsSymMatrix<D>>>;

  // Early return in case of range with length 1
  if (std::distance(begin, end) == 1) {
    return RetType{beginPars / beginWeight, *beginCov / beginWeight};
  }

  // Zero initialized values for aggregation
  ActsVector<D> mean = ActsVector<D>::Zero();
  ActsSymMatrix<D> cov1 = ActsSymMatrix<D>::Zero();
  ActsSymMatrix<D> cov2 = ActsSymMatrix<D>::Zero();
  WeightType sumOfWeights{0.0};
  bool haveCov = true;

  // clang-format off
  // x = \sum_{l} w_l * x_l
  // C = \sum_{l} w_l * C_l + \sum_{l} \sum_{m>l} w_l * w_m * (x_l - x_m)(x_l - x_m)^T
  // clang-format on
  for (auto l = begin; l != end; ++l) {
    const auto &[weight_l, pars_l, cov_l] = projector(*l);

    sumOfWeights += weight_l;
    mean += weight_l * pars_l;

    // Apply corrections for cyclic coordinates
    auto handleCyclicMean = [&ref = beginPars, &pars = pars_l,
                             &weight = weight_l, &mean = mean](auto desc) {
      const auto delta = (ref[desc.idx] - pars[desc.idx]) / desc.constant;

      if (delta > M_PI) {
        mean[desc.idx] += (2 * M_PI) * weight * desc.constant;
      } else if (delta < -M_PI) {
        mean[desc.idx] -= (2 * M_PI) * weight * desc.constant;
      }
    };

    std::apply([&](auto... dsc) { (handleCyclicMean(dsc), ...); }, angleDesc);

    if (not cov_l || not haveCov) {
      haveCov = false;
      continue;
    }

    cov1 += weight_l * *cov_l;

    for (auto m = std::next(l); m != end; ++m) {
      const auto &[weight_m, pars_m, cov_m] = projector(*m);

      ActsVector<D> diff = pars_l - pars_m;

      // Apply corrections for cyclic coordinates
      auto handleCyclicCov = [&l = pars_l, &m = pars_m,
                              &diff = diff](auto desc) {
        diff[desc.idx] =
            difference_periodic(l[desc.idx] / desc.constant,
                                m[desc.idx] / desc.constant, 2 * M_PI) *
            desc.constant;
      };

      std::apply([&](auto... dsc) { (handleCyclicCov(dsc), ...); }, angleDesc);

      cov2 += weight_l * weight_m * diff * diff.transpose();
    }
  }

  mean /= sumOfWeights;

  auto wrap = [&](auto desc) {
    mean[desc.idx] =
        wrap_periodic(mean[desc.idx] / desc.constant, -M_PI, 2 * M_PI) *
        desc.constant;
  };

  std::apply([&](auto... dsc) { (wrap(dsc), ...); }, angleDesc);

  if (haveCov) {
    cov1 /= sumOfWeights;
    cov2 /= sumOfWeights * sumOfWeights;
    return RetType{mean, cov1 + cov2};
  } else {
    return RetType{mean, std::nullopt};
  }
}

}  // namespace detail
}  // namespace Acts

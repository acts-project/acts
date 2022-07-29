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
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/Identity.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <cmath>
#include <optional>
#include <tuple>

namespace Acts {
namespace detail {

/// Description a cyclic angle
template <BoundIndices Idx>
struct CyclicAngle {
  constexpr static BoundIndices idx = Idx;
  constexpr static double constant = 1.0;
};

/// Description of a cyclic coordinate with a radius not known at compile time
template <BoundIndices Idx>
struct CyclicRadiusAngle {
  constexpr static BoundIndices idx = Idx;
  double constant = 1.0;  // the radius
};

/// A compile time map to provide angle descriptions for different surfaces
template <Surface::SurfaceType type_t>
struct AngleDescription {};

template <>
struct AngleDescription<Surface::Plane> {
  using Desc = std::tuple<CyclicAngle<eBoundPhi>>;
};

template <>
struct AngleDescription<Surface::Perigee> {
  using Desc = std::tuple<CyclicAngle<eBoundPhi>>;
};

template <>
struct AngleDescription<Surface::Disc> {
  using Desc = std::tuple<CyclicAngle<eBoundLoc1>, CyclicAngle<eBoundPhi>>;
};

template <>
struct AngleDescription<Surface::Cylinder> {
  using Desc =
      std::tuple<CyclicRadiusAngle<eBoundLoc0>, CyclicAngle<eBoundPhi>>;
};

/// A function object computing the approximate mean of a gaussian mixture
/// (approximations are applied for spherical, polar and angle coordinates)
struct CombineMeanApprox {
  template <typename components_t, typename projector_t, typename angle_desc_t>
  auto operator()(const components_t &cmps, const projector_t &proj,
                  const angle_desc_t &angleDesc) const {
    const auto &[beginWeight, beginPars, beginCov] = proj(cmps.front());

    constexpr int D = std::decay_t<decltype(beginPars)>::RowsAtCompileTime;

    ActsVector<D> mean = ActsVector<D>::Zero();
    double sumOfWeights = 0.0;

    for (const auto &cmp : cmps) {
      const auto &[weight, pars, cov] = proj(cmp);

      sumOfWeights += weight;
      mean += weight * pars;

      // Apply corrections for cyclic coordinates
      auto handleCyclicMean = [&ref = beginPars, &pars = pars, &weight = weight,
                               &mean = mean](auto desc) {
        const auto delta = (ref[desc.idx] - pars[desc.idx]) / desc.constant;

        if (delta > M_PI) {
          mean[desc.idx] += (2 * M_PI) * weight * desc.constant;
        } else if (delta < -M_PI) {
          mean[desc.idx] -= (2 * M_PI) * weight * desc.constant;
        }
      };

      std::apply([&](auto... dsc) { (handleCyclicMean(dsc), ...); }, angleDesc);
    }

    mean /= sumOfWeights;

    auto wrap = [&](auto desc) {
      mean[desc.idx] =
          wrap_periodic(mean[desc.idx] / desc.constant, -M_PI, 2 * M_PI) *
          desc.constant;
    };

    std::apply([&](auto... dsc) { (wrap(dsc), ...); }, angleDesc);

    return mean;
  }
};

/// Function object that computes the exact mean of Bound Coordinates on a plane
/// surface
template <Surface::SurfaceType type>
struct CombineMeanExact {
  template <typename components_t, typename projector_t, typename angle_desc_t>
  auto operator()(const components_t &cmps, const projector_t &proj,
                  const angle_desc_t &desc) const {
    // Basic mean for all cartesian coordinates
    BoundVector mean = CombineMeanApprox{}(cmps, proj, std::tuple<>{});

    // Direction mean in free coordinates
    Vector3 dir = Vector3::Zero();
    double sumOfWeights = 0.0;

    for (const auto &cmp : cmps) {
      const auto &[weight, pars, cov] = proj(cmp);

      sumOfWeights += weight;
      dir += weight *
             makeDirectionUnitFromPhiTheta(pars[eBoundPhi], pars[eBoundTheta]);
    }

    dir.normalize();
    mean.segment<2>(eBoundPhi) =
        Vector2{VectorHelpers::phi(dir), VectorHelpers::theta(dir)};

    // For plane surface we can return here
    if constexpr (type == Surface::Plane) {
      return mean;
    }

    // Special function for cylinder: circular mean
    auto circularMean = [&](auto i, auto c) {
      double x = 0.0, y = 0.0;

      for (const auto &cmp : cmps) {
        const auto &[weight, pars, cov] = proj(cmp);

        x += weight * std::cos(pars[i] / c);
        y += weight * std::sin(pars[i] / c);
      }

      return std::atan2(y / sumOfWeights, x / sumOfWeights) * c;
    };

    if constexpr (type == Surface::Cylinder) {
      mean[eBoundLoc0] = circularMean(eBoundLoc0, std::get<0>(desc).constant);
      return mean;
    }

    // Special function for disc and perigee: polar mean
    auto polarMean = [&](auto eR, auto ePhi) {
      double x = 0.0, y = 0.0;

      for (const auto &cmp : cmps) {
        const auto &[weight, pars, cov] = proj(cmp);
        x += weight * pars[eR] * std::cos(pars[ePhi]);
        y += weight * pars[eR] * std::sin(pars[ePhi]);
      }

      x /= sumOfWeights;
      y /= sumOfWeights;

      return Vector2{std::hypot(x, y), std::atan2(y, x)};
    };

    if constexpr (type == Surface::Disc) {
      mean.segment<2>(eBoundLoc0) = polarMean(eBoundLoc0, eBoundLoc1);
      return mean;
    }

    if constexpr (type == Surface::Perigee) {
      mean[eBoundLoc0] = polarMean(eBoundLoc0, eBoundPhi)[0];
      return mean;
    }
  }
};

template <int D, typename components_t, typename projector_t,
          typename angle_desc_t>
auto combineCov(const components_t &components, const ActsVector<D> &mean,
                const projector_t &proj, const angle_desc_t &angleDesc) {
  ActsSymMatrix<D> aggregatedCov = ActsSymMatrix<D>::Zero();
  double sumOfWeights{0.0};

  for (const auto &cmp : components) {
    const auto &[weight, pars, cov] = proj(cmp);

    sumOfWeights += weight;
    ActsVector<D> diff = pars - mean;

    // Apply corrections for cyclic coordinates
    auto handleCyclicCov = [&l = pars, &m = mean, &diff = diff](auto desc) {
      diff[desc.idx] =
          difference_periodic(l[desc.idx] / desc.constant,
                              m[desc.idx] / desc.constant, 2 * M_PI) *
          desc.constant;
    };

    std::apply([&](auto... dsc) { (handleCyclicCov(dsc), ...); }, angleDesc);

    aggregatedCov += weight * (*cov + diff * diff.transpose());
  }

  aggregatedCov /= sumOfWeights;
  return aggregatedCov;
}

/// @brief Combine multiple components into one representative track state
/// object. The function takes iterators to allow for arbitrary ranges to be
/// combined. The dimension of the vectors is infeared from the inputs.
///
/// @note The behaviour is undefined, if some components have a covariance
/// and others not
///
/// @note The correct mean and variances for cyclic coordnates or spherical
/// coordinates (theta, phi) must generally be computed using a special circular
/// mean or in cartesian coordinates. This implements a approximation, which
/// only works well for close components.
///
/// @tparam component_iterator_t An iterator of a range of components
/// @tparam projector_t A proj, which maps the component to a
/// std::tuple< weight, mean, std::optional< cov > >
/// @tparam angle_desc_t A angle description object which defines the cyclic
/// angles in the bound parameters
template <typename components_t, typename projector_t = Identity,
          typename angle_desc_t = AngleDescription<Surface::Plane>::Desc,
          typename combine_mean_t = CombineMeanApprox>
auto combineGaussianMixture(
    const components_t &components, projector_t &&proj = projector_t{},
    const angle_desc_t &angleDesc = angle_desc_t{},
    const combine_mean_t &combineMean = combine_mean_t{}) {
  // Extract the first component
  const auto &[beginWeight, beginPars, beginCov] = proj(components.front());

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
  if (components.size() == 1) {
    return RetType{beginPars / beginWeight, *beginCov / beginWeight};
  }

  const auto mean = combineMean(components, proj, angleDesc);

  if (beginCov) {
    const auto cov = combineCov(components, mean, proj, angleDesc);
    return RetType{std::move(mean), std::move(cov)};
  } else {
    return RetType{std::move(mean), std::nullopt};
  }
}

template <bool exact, Surface::SurfaceType type>
auto makeCombineFunction(const Surface &surface) {
  // Make angle description object
  auto desc = typename AngleDescription<type>::Desc{};

  if constexpr (type == Surface::Cylinder) {
    const auto &bounds = static_cast<const CylinderSurface &>(surface).bounds();
    std::get<0>(desc).constant = bounds.get(CylinderBounds::eR);
  }

  // Exact functions
  if constexpr (exact) {
    return [=](const auto &cmps, auto &&proj) {
      return combineGaussianMixture(cmps, proj, desc, CombineMeanExact<type>{});
    };
  } else {
    return [=](const auto &cmps, auto &&proj) {
      return combineGaussianMixture(cmps, proj, desc, CombineMeanApprox{});
    };
  }
}

}  // namespace detail
}  // namespace Acts

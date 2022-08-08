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
#include <numeric>
#include <optional>
#include <tuple>

#include <Eigen/Eigenvalues>

namespace Acts {
namespace detail {

/// Angle descriptions for the combineBoundGaussianMixture function
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

/// A compile time map to provide angle descriptions for different surfaces
template <Surface::SurfaceType type_t>
struct AngleDescription {
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

// Helper function that encapsulates a switch-case to select the correct angle
// description dependent on the surface
template <typename Callable>
auto angleDescriptionSwitch(const Surface &surface, Callable &&callable) {
  switch (surface.type()) {
    case Surface::Cylinder: {
      auto desc = AngleDescription<Surface::Cylinder>::Desc{};
      const auto &bounds =
          static_cast<const CylinderSurface &>(surface).bounds();
      std::get<0>(desc).constant = bounds.get(CylinderBounds::eR);
      return callable(desc);
    }
    case Surface::Disc: {
      auto desc = AngleDescription<Surface::Disc>::Desc{};
      return callable(desc);
    }
    default: {
      auto desc = AngleDescription<Surface::Plane>::Desc{};
      return callable(desc);
    }
  }
}

template <int D, typename components_t, typename projector_t,
          typename angle_desc_t>
auto combineCov(const components_t components, const ActsVector<D> &mean,
                double sumOfWeights, projector_t &&projector,
                const angle_desc_t &angleDesc) {
  ActsSymMatrix<D> cov = ActsSymMatrix<D>::Zero();

  for (const auto &cmp : components) {
    const auto &[weight_l, pars_l, cov_l] = projector(cmp);

    cov += weight_l * *cov_l;

    ActsVector<D> diff = pars_l - mean;

    // Apply corrections for cyclic coordinates
    auto handleCyclicCov = [&l = pars_l, &m = mean, &diff = diff](auto desc) {
      diff[desc.idx] =
          difference_periodic(l[desc.idx] / desc.constant,
                              m[desc.idx] / desc.constant, 2 * M_PI) *
          desc.constant;
    };

    std::apply([&](auto... dsc) { (handleCyclicCov(dsc), ...); }, angleDesc);

    cov += weight_l * diff * diff.transpose();
  }

  cov /= sumOfWeights;
  return cov;
}

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
/// @tparam components_t A range of components
/// @tparam projector_t A projector, which maps the component to a
/// std::tuple< weight, mean, std::optional< cov > >
/// @tparam angle_desc_t A angle description object which defines the cyclic
/// angles in the bound parameters
template <typename components_t, typename projector_t = Identity,
          typename angle_desc_t = AngleDescription<Surface::Plane>::Desc>
auto combineGaussianMixture(const components_t components,
                            projector_t &&projector = projector_t{},
                            const angle_desc_t &angleDesc = angle_desc_t{}) {
  // Extract the first component
  const auto &[beginWeight, beginPars, beginCov] =
      projector(components.front());

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

  // Zero initialized values for aggregation
  ActsVector<D> mean = ActsVector<D>::Zero();
  WeightType sumOfWeights{0.0};

  for (const auto &cmp : components) {
    const auto &[weight_l, pars_l, cov_l] = projector(cmp);

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
  }

  mean /= sumOfWeights;

  auto wrap = [&](auto desc) {
    mean[desc.idx] =
        wrap_periodic(mean[desc.idx] / desc.constant, -M_PI, 2 * M_PI) *
        desc.constant;
  };

  std::apply([&](auto... dsc) { (wrap(dsc), ...); }, angleDesc);

  if (not beginCov) {
    return RetType{mean, std::nullopt};
  } else {
    const auto cov =
        combineCov(components, mean, sumOfWeights, projector, angleDesc);

    return RetType{mean, cov};
  }
}


template <typename components_t, typename projector_t>
auto unevaluated_context_mean(const components_t &cmps, const projector_t &proj)
{
  const auto &[w, mean, c] = proj(cmps.front());
  return mean;
}

/// This function implements a mode-search using an iterative fixed-point algorithm
/// TODO potential optimization: cache inverse covariances
// https://faculty.ucmerced.edu/mcarreira-perpinan/papers/cs-99-03.pdf
template <typename components_t, typename projector_t = Identity,
          typename angle_desc_t = AngleDescription<Surface::Plane>::Desc>
auto computeModeOfMixture(const components_t components,
                            projector_t &&proj = projector_t{},
                            const angle_desc_t &desc = angle_desc_t{}) {
  constexpr int D = [](){
    using ParsType = decltype(unevaluated_context_mean(std::declval<components_t>(), std::declval<projector_t>()));
    return ParsType::RowsAtCompileTime;
  }();

  // Lambda used to correct cyclic coordinates in the computation
  auto cyclicDiff = [](auto d, auto &diff, const auto &a, const auto &b) {
    diff[d.idx] =
        difference_periodic(a[d.idx] / d.constant,
                            b[d.idx] / d.constant, 2 * M_PI) *
        d.constant;
  };

  // Compute the value of the pdf of a single multivariate gaussian
  auto single_pdf = [&](const auto &x, const auto &mean, const auto &cov) {
    const auto a = 1.0 / std::sqrt(std::pow(2 * M_PI, D) * cov->determinant());
    ActsVector<D> r = x - mean;
    std::apply([&](auto... d) { (cyclicDiff(d, r, x, mean), ...); }, desc);
    const auto b = -0.5 * (x - mean).transpose() * cov->inverse() * (x - mean);
    return a * std::exp(b);
  };

  // Compute the value of the gaussian mixture pdf at x
  auto mixture_pdf = [&](const auto &x) {
    double res = 0.0;

    for (const auto &cmp : components) {
      const auto &[weight, mean, cov] = proj(cmp);
      res += weight * single_pdf(x, mean, cov);
    }

    return res;
  };

  // The fixed-point equation (see ref)
  auto f = [&](const auto &x) {
    const auto p_x = mixture_pdf(x);

    ActsSymMatrix<D> a = ActsSymMatrix<D>::Zero();
    ActsVector<D> b = ActsVector<D>::Zero();

    for (const auto &cmp : components) {
      const auto &[weight, mean, cov] = proj(cmp);
      const auto p_m_x = weight * single_pdf(x, mean, cov) / p_x;

      a += p_m_x * cov->inverse();
      b += p_m_x * (cov->inverse() * mean);
    }

    return a.inverse() * b;
  };

  // Compute the hessian of a gaussian mixture at x
  auto hessian = [&](const auto &x) {
    ActsSymMatrix<D> h = ActsSymMatrix<D>::Zero();

    for (const auto &cmp : components) {
      const auto &[weight, mean, cov] = proj(cmp);
      const auto a = weight * single_pdf(x, mean, cov);
      ActsVector<D> r = x - mean;
      std::apply([&](auto... d) { (cyclicDiff(d, r, x, mean), ...); },
                 desc);
      const auto b = (x - mean) * (x - mean).transpose() - *cov;
      h += a * (cov->inverse() * b * cov->inverse());
    }

    return h;
  };

  // We only store the highest mode
  std::optional<ActsVector<D>> mode;
  double mode_val = 0;

  // Algorithm parameters
  constexpr double tol = 1.e-8;
  constexpr double eigv_max = 0.01;
  constexpr int iter_max = 50;

  // The algorithm
  for (const auto &cmp : components) {
    const auto &[weight, mean, cov] = proj(cmp);

    ActsVector<D> x = mean;
    ActsVector<D> x_old = ActsVector<D>::Zero();

    for (auto i = 0; i < iter_max; ++i) {
      x_old = x;
      x = f(x);

      if ((x - x_old).norm() < tol) {
        break;
      }
    }

    const auto H = hessian(x);
    const auto evs = H.eigenvalues();
    const double max_ev = std::max_element(evs.data(), evs.data() + evs.size(),
                                           [](const auto &a, const auto &b) {
                                             return a.real() < b.real();
                                           })
                              ->real();

    if (max_ev < eigv_max) {
      const auto this_val = mixture_pdf(x);
      if (this_val > mode_val) {
        mode = x;
        mode_val = this_val;
      }
    }
  }

  return mode;
}

}  // namespace detail
}  // namespace Acts

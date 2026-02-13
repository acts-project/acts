// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/TrackFitting/GsfComponent.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <iosfwd>
#include <stdexcept>
#include <tuple>

namespace Acts::detail::Gsf {

/// Reduce Gaussian mixture
/// @param mixture The mixture iterable
/// @param surface The surface, on which the bound state is
/// @param method How to reduce the mixture
/// @return parameters and covariance as tuple
std::tuple<BoundVector, BoundMatrix> mergeGaussianMixture(
    std::span<const GsfComponent> mixture, const Surface &surface,
    ComponentMergeMethod method);

/// Merge two GSF components into one
/// @param a First component
/// @param b Second component
/// @param surface The surface, on which the bound state is
/// @return merged component
GsfComponent mergeTwoComponents(const GsfComponent &a, const GsfComponent &b,
                                const Surface &surface);

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

/// Helper function that encapsulates a switch-case to select the correct angle
/// description dependent on the surface
/// @param surface The surface, which the bound state is on
/// @param callable A callable that takes the angle description as argument
/// @return The return value of the callable
template <typename Callable>
auto angleDescriptionSwitch(const Surface &surface, Callable &&callable) {
  switch (surface.type()) {
    case Surface::Cylinder: {
      AngleDescription<Surface::Cylinder>::Desc desc{};
      const auto &bounds =
          static_cast<const CylinderSurface &>(surface).bounds();
      std::get<0>(desc).constant = bounds.get(CylinderBounds::eR);
      return callable(desc);
    }
    case Surface::Disc: {
      AngleDescription<Surface::Disc>::Desc desc{};
      return callable(desc);
    }
    default: {
      AngleDescription<Surface::Plane>::Desc desc{};
      return callable(desc);
    }
  }
}

/// Compute the covariance of a Gaussian mixture given the mean. The function
/// takes iterators to allow for arbitrary ranges to be combined. The dimension
/// of the vectors is inferred from the inputs.
///
/// @param cmps The component range to merge
/// @param projector A projector to extract the weight, parameters and covariance
///        from the components
/// @param mean The precomputed mean of the mixture
/// @param sumOfWeights The precomputed sum of weights of the mixture
/// @param angleDesc The angle description object
/// @return covariance
template <typename component_range_t, typename projector_t,
          typename angle_desc_t>
BoundMatrix mergeGaussianMixtureCov(const component_range_t &cmps,
                                    const projector_t &projector,
                                    const BoundVector &mean,
                                    double sumOfWeights,
                                    const angle_desc_t &angleDesc) {
  BoundMatrix cov = BoundMatrix::Zero();

  for (const auto &cmp : cmps) {
    const auto &[weight_l, pars_l, cov_l] = projector(cmp);

    cov += weight_l * cov_l;  // MARK: fpeMask(FLTUND, 1, #2347)

    BoundVector diff = pars_l - mean;

    // Apply corrections for cyclic coordinates
    const auto handleCyclicCov = [&l = pars_l, &m = mean,
                                  &diff = diff](const auto &desc) {
      diff[desc.idx] = difference_periodic(l[desc.idx] / desc.constant,
                                           m[desc.idx] / desc.constant,
                                           2 * std::numbers::pi) *
                       desc.constant;
    };
    std::apply([&](auto... dsc) { (handleCyclicCov(dsc), ...); }, angleDesc);

    cov += weight_l * diff * diff.transpose();
  }

  cov /= sumOfWeights;
  return cov;
}

/// Combine multiple components into one representative track state object. The
/// function takes iterators to allow for arbitrary ranges to be combined. The
/// dimension of the vectors is inferred from the inputs.
///
/// @note If one component does not contain a covariance, no covariance is
/// computed.
///
/// @note The correct mean and variances for cyclic coordinates or spherical
/// coordinates (theta, phi) must generally be computed using a special circular
/// mean or in cartesian coordinates. This implements a approximation, which
/// only works well for close components.
///
/// @param cmps The component range to merge
/// @param projector A projector to extract the weight, parameters and covariance
///        from the components
/// @param angleDesc The angle description object
/// @return parameters and covariance
template <typename component_range_t, typename projector_t,
          typename angle_desc_t>
std::tuple<BoundVector, BoundMatrix> mergeGaussianMixtureMeanCov(
    const component_range_t &cmps, const projector_t &projector,
    const angle_desc_t &angleDesc) {
  // Early return in case of range with length 1
  if (cmps.size() == 1) {
    const auto &[beginWeight, beginPars, beginCov] = projector(cmps.front());
    return {beginPars / beginWeight, beginCov / beginWeight};
  }

  // Do the (circular) mean with complex arithmetic.
  // For normal values, just keep the real values. For angles, use the complex
  // phase. Weighting then transparently happens by multiplying a real-valued
  // weight.
  using CVec = Eigen::Matrix<std::complex<double>, eBoundSize, 1>;
  CVec cMean = CVec::Zero();
  double sumOfWeights{0.0};

  for (const auto &cmp : cmps) {
    const auto &[weight_l, pars_l, cov_l] = projector(cmp);

    CVec pars_l_c = pars_l;

    const auto setPolar = [&](const auto &desc) {
      pars_l_c[desc.idx] = std::polar(1.0, pars_l[desc.idx] / desc.constant);
    };
    std::apply([&](auto... dsc) { (setPolar(dsc), ...); }, angleDesc);

    sumOfWeights += weight_l;
    cMean += weight_l * pars_l_c;
  }

  cMean /= sumOfWeights;

  BoundVector mean = cMean.real();

  const auto getArg = [&](const auto &desc) {
    mean[desc.idx] = desc.constant * std::arg(cMean[desc.idx]);
  };
  std::apply([&](auto... dsc) { (getArg(dsc), ...); }, angleDesc);

  // MARK: fpeMaskBegin(FLTUND, 1, #2347)
  const BoundMatrix cov =
      mergeGaussianMixtureCov(cmps, projector, mean, sumOfWeights, angleDesc);
  // MARK: fpeMaskEnd(FLTUND)

  return {mean, cov};
}

/// Reduce Gaussian mixture
///
/// @param cmps The component range to merge
/// @param projector A projector to extract the weight, parameters and covariance
///        from the components
/// @param angleDesc The angle description object
/// @param method How to reduce the mixture
/// @return parameters and covariance
template <typename component_range_t, typename projector_t,
          typename angle_desc_t>
std::tuple<BoundVector, BoundMatrix> mergeGaussianMixture(
    const component_range_t &cmps, const projector_t &projector,
    const angle_desc_t &angleDesc, ComponentMergeMethod method) {
  const auto [mean, cov] =
      mergeGaussianMixtureMeanCov(cmps, projector, angleDesc);

  if (method == ComponentMergeMethod::eMean) {
    return {mean, cov};
  } else if (method == ComponentMergeMethod::eMaxWeight) {
    const auto maxWeightIt =
        std::ranges::max_element(cmps, {}, [&](const auto &cmp) {
          const auto &[weight_l, pars_l, cov_l] = projector(cmp);
          return weight_l;
        });
    const auto &[weight_l, pars_l, cov_l] = projector(*maxWeightIt);

    return {pars_l, cov};
  } else {
    throw std::logic_error("Invalid component merge method");
  }
}

/// Reduce Gaussian mixture
///
/// @param cmps The component range to merge
/// @param projector A projector to extract the weight, parameters and covariance
///        from the components
/// @param surface The surface, which the bound state is on
/// @param method How to reduce the mixture
/// @return parameters and covariance
template <typename component_range_t, typename projector_t>
std::tuple<BoundVector, BoundMatrix> mergeGaussianMixture(
    const component_range_t &cmps, const projector_t &projector,
    const Surface &surface, ComponentMergeMethod method) {
  return angleDescriptionSwitch(surface, [&](const auto &desc) {
    return mergeGaussianMixture(cmps, projector, desc, method);
  });
}

/// @brief Class representing a symmetric distance matrix
class SymmetricKLDistanceMatrix {
  using Array = Eigen::Array<double, Eigen::Dynamic, 1>;
  using Mask = Eigen::Array<bool, Eigen::Dynamic, 1>;

  Array m_distances;
  Mask m_mask;
  std::vector<std::pair<std::size_t, std::size_t>> m_mapToPair;
  std::size_t m_numberComponents{};

  template <typename array_t, typename setter_t>
  void setAssociated(std::size_t n, array_t &array, const setter_t &setter) {
    const std::size_t indexConst = (n - 1) * n / 2;
    // Rows
    for (std::size_t i = 0ul; i < n; ++i) {
      array[indexConst + i] = setter(n, i);
    }
    // Columns
    for (std::size_t i = n + 1; i < m_numberComponents; ++i) {
      array[(i - 1) * i / 2 + n] = setter(n, i);
    }
  }

  /// Computes the Kullback-Leibler distance between two components as shown in
  /// https://arxiv.org/abs/2001.00727v1 but ignoring the weights
  static double computeSymmetricKlDivergence(const GsfComponent &a,
                                             const GsfComponent &b);

 public:
  explicit SymmetricKLDistanceMatrix(std::span<const GsfComponent> cmps);

  double at(std::size_t i, std::size_t j) const;

  void recomputeAssociatedDistances(std::size_t n,
                                    std::span<const GsfComponent> cmps);

  void maskAssociatedDistances(std::size_t n);

  std::pair<std::size_t, std::size_t> minDistancePair() const;

  std::ostream &toStream(std::ostream &os) const;

  friend std::ostream &operator<<(std::ostream &os,
                                  const SymmetricKLDistanceMatrix &m) {
    return m.toStream(os);
  }
};

}  // namespace Acts::detail::Gsf

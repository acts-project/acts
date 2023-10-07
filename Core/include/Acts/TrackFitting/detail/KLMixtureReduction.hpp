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
#include "Acts/Utilities/detail/gaussian_mixture_helpers.hpp"

#include "GsfUtils.hpp"

namespace Acts {

namespace detail {

/// Computes the Kullback-Leibler distance between two components as shown in
/// https://arxiv.org/abs/2001.00727v1 but ignoring the weights
template <typename component_t, typename component_projector_t>
auto computeSymmetricKlDivergence(const component_t &a, const component_t &b,
                                  const component_projector_t &proj) {
  const auto parsA = proj(a).boundPars[eBoundQOverP];
  const auto parsB = proj(b).boundPars[eBoundQOverP];
  const auto covA = proj(a).boundCov(eBoundQOverP, eBoundQOverP);
  const auto covB = proj(b).boundCov(eBoundQOverP, eBoundQOverP);

  throw_assert(covA != 0.0, "");
  throw_assert(std::isfinite(covA), "");
  throw_assert(covB != 0.0, "");
  throw_assert(std::isfinite(covB), "");

  const auto kl = covA * (1 / covB) + covB * (1 / covA) +
                  (parsA - parsB) * (1 / covA + 1 / covB) * (parsA - parsB);

  throw_assert(kl >= 0.0, "kl-distance should be positive, but is: "
                              << kl << "(qop_a: " << parsA << "+-" << covA
                              << ", qop_b: " << parsB << "+-" << covB << ")");

  return kl;
}

template <typename component_t, typename component_projector_t,
          typename angle_desc_t>
auto mergeComponents(const component_t &a, const component_t &b,
                     const component_projector_t &proj,
                     const angle_desc_t &angle_desc) {
  throw_assert(proj(a).weight > 0.0 && proj(b).weight > 0.0, "weight error");

  std::array range = {std::ref(proj(a)), std::ref(proj(b))};
  const auto refProj = [](auto &c) {
    return std::tie(c.get().weight, c.get().boundPars, c.get().boundCov);
  };

  auto [mergedPars, mergedCov] =
      combineGaussianMixture(range, refProj, angle_desc);

  component_t ret = a;
  proj(ret).boundPars = mergedPars;
  proj(ret).boundCov = mergedCov;
  proj(ret).weight = proj(a).weight + proj(b).weight;

  return ret;
}

/// @brief Class representing a symmetric distance matrix
class SymmetricKLDistanceMatrix {
  using Array = Eigen::Array<ActsScalar, Eigen::Dynamic, 1>;
  using Mask = Eigen::Array<bool, Eigen::Dynamic, 1>;

  Array m_distances;
  Mask m_mask;
  std::vector<std::pair<std::size_t, std::size_t>> m_mapToPair;
  std::size_t m_numberComponents;

  template <typename array_t, typename setter_t>
  void setAssociated(std::size_t n, array_t &array, setter_t &&setter) {
    const auto indexConst = (n - 1) * n / 2;

    // Rows
    for (auto i = 0ul; i < n; ++i) {
      array[indexConst + i] = setter(n, i);
    }

    // Columns
    for (auto i = n + 1; i < m_numberComponents; ++i) {
      array[(i - 1) * i / 2 + n] = setter(n, i);
    }
  }

 public:
  template <typename component_t, typename projector_t>
  SymmetricKLDistanceMatrix(const std::vector<component_t> &cmps,
                            const projector_t &proj)
      : m_distances(Array::Zero(cmps.size() * (cmps.size() - 1) / 2)),
        m_mask(Mask::Ones(cmps.size() * (cmps.size() - 1) / 2)),
        m_mapToPair(m_distances.size()),
        m_numberComponents(cmps.size()) {
    for (auto i = 1ul; i < m_numberComponents; ++i) {
      const auto indexConst = (i - 1) * i / 2;
      for (auto j = 0ul; j < i; ++j) {
        m_mapToPair.at(indexConst + j) = {i, j};
        m_distances[indexConst + j] =
            computeSymmetricKlDivergence(cmps[i], cmps[j], proj);
      }
    }
  }

  auto at(std::size_t i, std::size_t j) const {
    return m_distances[i * (i - 1) / 2 + j];
  }

  template <typename component_t, typename projector_t>
  void recomputeAssociatedDistances(std::size_t n,
                                    const std::vector<component_t> &cmps,
                                    const projector_t &proj) {
    assert(cmps.size() == m_numberComponents && "size mismatch");

    setAssociated(n, m_distances, [&](std::size_t i, std::size_t j) {
      return computeSymmetricKlDivergence(cmps[i], cmps[j], proj);
    });
  }

  void maskAssociatedDistances(std::size_t n) {
    setAssociated(n, m_mask, [&](std::size_t, std::size_t) { return false; });
  }

  auto minDistancePair() const {
    ActsScalar min = std::numeric_limits<ActsScalar>::max();
    std::size_t idx = 0;

    for (auto i = 0l; i < m_distances.size(); ++i) {
      if (auto new_min = std::min(min, m_distances[i]);
          m_mask[i] && new_min < min) {
        min = new_min;
        idx = i;
      }
    }

    return m_mapToPair.at(idx);
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const SymmetricKLDistanceMatrix &m) {
    const auto prev_precision = os.precision();
    const int width = 8;
    const int prec = 2;

    os << "\n";
    os << std::string(width, ' ') << " | ";
    for (auto j = 0ul; j < m.m_numberComponents - 1; ++j) {
      os << std::setw(width) << j << "  ";
    }
    os << "\n";
    os << std::string((width + 3) + (width + 2) * (m.m_numberComponents - 1),
                      '-');
    os << "\n";

    for (auto i = 1ul; i < m.m_numberComponents; ++i) {
      const auto indexConst = (i - 1) * i / 2;
      os << std::setw(width) << i << " | ";
      for (auto j = 0ul; j < i; ++j) {
        os << std::setw(width) << std::setprecision(prec)
           << m.m_distances[indexConst + j] << "  ";
      }
      os << "\n";
    }
    os << std::setprecision(prev_precision);
    return os;
  }
};

template <typename component_t, typename component_projector_t,
          typename angle_desc_t>
void reduceWithKLDistance(std::vector<component_t> &cmpCache,
                          std::size_t maxCmpsAfterMerge,
                          const component_projector_t &proj,
                          const angle_desc_t &angle_desc) {
  if (cmpCache.size() <= maxCmpsAfterMerge) {
    return;
  }

  SymmetricKLDistanceMatrix distances(cmpCache, proj);

  auto remainingComponents = cmpCache.size();

  while (remainingComponents > maxCmpsAfterMerge) {
    const auto [minI, minJ] = distances.minDistancePair();

    // Set one component and compute associated distances
    cmpCache[minI] =
        mergeComponents(cmpCache[minI], cmpCache[minJ], proj, angle_desc);
    distances.recomputeAssociatedDistances(minI, cmpCache, proj);

    // Set weight of the other component to -1 so we can remove it later and
    // mask its distances
    proj(cmpCache[minJ]).weight = -1.0;
    distances.maskAssociatedDistances(minJ);

    remainingComponents--;
  }

  // Remove all components which are labled with weight -1
  std::sort(cmpCache.begin(), cmpCache.end(),
            [&](const auto &a, const auto &b) {
              return proj(a).weight < proj(b).weight;
            });
  cmpCache.erase(
      std::remove_if(cmpCache.begin(), cmpCache.end(),
                     [&](const auto &a) { return proj(a).weight == -1.0; }),
      cmpCache.end());

  throw_assert(cmpCache.size() == maxCmpsAfterMerge,
               "size mismatch, should be " << maxCmpsAfterMerge << ", but is "
                                           << cmpCache.size());
}

}  // namespace detail

}  // namespace Acts

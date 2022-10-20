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
/// https://arxiv.org/abs/2001.00727v1, while ignoring the weights as done in
/// Athena
template <typename component_t, typename component_projector_t>
auto computeKLDistance(const component_t &a, const component_t &b,
                       const component_projector_t &proj) {
  const auto parsA = proj(a).boundPars[eBoundQOverP];
  const auto parsB = proj(b).boundPars[eBoundQOverP];
  const auto covA = (*proj(a).boundCov)(eBoundQOverP, eBoundQOverP);
  const auto covB = (*proj(b).boundCov)(eBoundQOverP, eBoundQOverP);

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
  Eigen::VectorXd m_data;
  std::vector<std::pair<std::size_t, std::size_t>> m_mapToPair;
  std::size_t m_N;

 public:
  template <typename component_t, typename projector_t>
  SymmetricKLDistanceMatrix(const std::vector<component_t> &cmps,
                            const projector_t &proj)
      : m_data(cmps.size() * (cmps.size() - 1) / 2),
        m_mapToPair(m_data.size()),
        m_N(cmps.size()) {
    for (auto i = 1ul; i < m_N; ++i) {
      const auto indexConst = (i - 1) * i / 2;
      for (auto j = 0ul; j < i; ++j) {
        m_mapToPair.at(indexConst + j) = {i, j};
        m_data[indexConst + j] = computeKLDistance(cmps[i], cmps[j], proj);
      }
    }
  }

  auto at(std::size_t i, std::size_t j) const {
    return m_data[i * (i - 1) / 2 + j];
  }

  template <typename component_t, typename projector_t>
  void recomputeAssociatedDistances(std::size_t n,
                                    const std::vector<component_t> &cmps,
                                    const projector_t &proj) {
    const auto indexConst = (n - 1) * n / 2;

    throw_assert(cmps.size() == m_N, "size mismatch");

    // Rows
    for (auto i = 0ul; i < n; ++i) {
      m_data[indexConst + i] = computeKLDistance(cmps[n], cmps[i], proj);
    }

    // Columns
    for (auto i = n + 1; i < cmps.size(); ++i) {
      m_data[(i - 1) * i / 2 + n] = computeKLDistance(cmps[n], cmps[i], proj);
    }
  }

  void resetAssociatedDistances(std::size_t n, double value) {
    const auto indexConst = (n - 1) * n / 2;

    // Rows
    for (auto i = 0ul; i < n; ++i) {
      m_data[indexConst + i] = value;
    }

    // Columns
    for (auto i = n + 1; i < m_N; ++i) {
      m_data[(i - 1) * i / 2 + n] = value;
    }
  }

  auto minDistancePair() const {
    // TODO Eigen minCoeff does not work for some reason??
    // return m_mapToPair.at(m_data.minCoeff());
    return m_mapToPair.at(std::distance(
        &m_data[0],
        std::min_element(m_data.data(), m_data.data() + m_data.size())));
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const SymmetricKLDistanceMatrix &m) {
    for (auto i = 1ul; i < m.m_N; ++i) {
      const auto indexConst = (i - 1) * i / 2;
      Eigen::RowVectorXd vals;
      vals.resize(i);
      for (auto j = 0ul; j < i; ++j) {
        vals[j] = m.m_data[indexConst + j];
      }
      os << vals << "\n";
    }

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

    cmpCache[minI] =
        mergeComponents(cmpCache[minI], cmpCache[minJ], proj, angle_desc);
    distances.recomputeAssociatedDistances(minI, cmpCache, proj);
    remainingComponents--;

    // Reset removed components so that it won't have the shortest distance
    // ever, and so that we can sort them by weight in the end to remove them
    proj(cmpCache[minJ]).weight = -1.0;
    proj(cmpCache[minJ]).boundPars[eBoundQOverP] =
        std::numeric_limits<double>::max();
    (*proj(cmpCache[minJ]).boundCov)(eBoundQOverP, eBoundQOverP) =
        std::numeric_limits<double>::max();
    distances.resetAssociatedDistances(minJ,
                                       std::numeric_limits<double>::max());
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

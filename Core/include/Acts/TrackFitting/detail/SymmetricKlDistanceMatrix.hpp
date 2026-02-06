// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFitting/GsfComponent.hpp"

#include <iomanip>
#include <span>

namespace Acts::detail {

/// Computes the Kullback-Leibler distance between two components as shown in
/// https://arxiv.org/abs/2001.00727v1 but ignoring the weights
auto computeSymmetricKlDivergence(const GsfComponent &a,
                                  const GsfComponent &b) {
  const auto parsA = a.boundPars[eBoundQOverP];
  const auto parsB = b.boundPars[eBoundQOverP];
  const auto covA = a.boundCov(eBoundQOverP, eBoundQOverP);
  const auto covB = b.boundCov(eBoundQOverP, eBoundQOverP);

  assert(covA != 0.0);
  assert(std::isfinite(covA));
  assert(covB != 0.0);
  assert(std::isfinite(covB));

  const auto kl = covA * (1 / covB) + covB * (1 / covA) +
                  (parsA - parsB) * (1 / covA + 1 / covB) * (parsA - parsB);

  assert(kl >= 0.0 && "kl-divergence must be non-negative");

  return kl;
}

/// @brief Class representing a symmetric distance matrix
class SymmetricKLDistanceMatrix {
  using Array = Eigen::Array<double, Eigen::Dynamic, 1>;
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
  explicit SymmetricKLDistanceMatrix(std::span<const GsfComponent> cmps)
      : m_distances(Array::Zero(cmps.size() * (cmps.size() - 1) / 2)),
        m_mask(Mask::Ones(cmps.size() * (cmps.size() - 1) / 2)),
        m_mapToPair(m_distances.size()),
        m_numberComponents(cmps.size()) {
    for (auto i = 1ul; i < m_numberComponents; ++i) {
      const auto indexConst = (i - 1) * i / 2;
      for (auto j = 0ul; j < i; ++j) {
        m_mapToPair.at(indexConst + j) = {i, j};
        m_distances[indexConst + j] =
            computeSymmetricKlDivergence(cmps[i], cmps[j]);
      }
    }
  }

  auto at(std::size_t i, std::size_t j) const {
    return m_distances[i * (i - 1) / 2 + j];
  }

  void recomputeAssociatedDistances(std::size_t n,
                                    std::span<const GsfComponent> cmps) {
    assert(cmps.size() == m_numberComponents && "size mismatch");

    setAssociated(n, m_distances, [&](std::size_t i, std::size_t j) {
      return computeSymmetricKlDivergence(cmps[i], cmps[j]);
    });
  }

  void maskAssociatedDistances(std::size_t n) {
    setAssociated(n, m_mask, [&](std::size_t, std::size_t) { return false; });
  }

  auto minDistancePair() const {
    auto min = std::numeric_limits<double>::max();
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

}  // namespace Acts::detail

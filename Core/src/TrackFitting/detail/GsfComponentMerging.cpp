// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/detail/GsfComponentMerging.hpp"

#include <algorithm>
#include <format>
#include <iostream>
#include <limits>
#include <numeric>

namespace Acts {

std::tuple<BoundVector, BoundMatrix> detail::Gsf::mergeGaussianMixture(
    std::span<const GsfComponent> mixture, const Surface &surface,
    ComponentMergeMethod method) {
  return mergeGaussianMixture(
      mixture,
      [](const GsfComponent &c) {
        return std::tie(c.weight, c.boundPars, c.boundCov);
      },
      surface, method);
}

GsfComponent detail::Gsf::mergeTwoComponents(const GsfComponent &a,
                                             const GsfComponent &b,
                                             const Surface &surface) {
  assert(a.weight >= 0.0 && b.weight >= 0.0 && "non-positive weight");

  std::array components = {&a, &b};
  const auto proj = [](const GsfComponent *c) {
    return std::tie(c->weight, c->boundPars, c->boundCov);
  };
  auto [mergedPars, mergedCov] =
      angleDescriptionSwitch(surface, [&](const auto &desc) {
        return mergeGaussianMixtureMeanCov(components, proj, desc);
      });

  GsfComponent ret = a;
  ret.boundPars = mergedPars;
  ret.boundCov = mergedCov;
  ret.weight = a.weight + b.weight;
  return ret;
}

double detail::Gsf::computeSymmetricKlDivergence(const GsfComponent &a,
                                                 const GsfComponent &b) {
  const double parsA = a.boundPars[eBoundQOverP];
  const double parsB = b.boundPars[eBoundQOverP];
  const double covA = a.boundCov(eBoundQOverP, eBoundQOverP);
  const double covB = b.boundCov(eBoundQOverP, eBoundQOverP);

  assert(covA != 0.0);
  assert(std::isfinite(covA));
  assert(covB != 0.0);
  assert(std::isfinite(covB));

  const double kl = covA * (1 / covB) + covB * (1 / covA) +
                    (parsA - parsB) * (1 / covA + 1 / covB) * (parsA - parsB);

  assert(kl >= 0.0 && "kl-divergence must be non-negative");

  return kl;
}

namespace detail::Gsf {

SymmetricKLDistanceMatrix::SymmetricKLDistanceMatrix(
    std::span<const GsfComponent> cmps)
    : m_distances(Array::Zero(cmps.size() * (cmps.size() - 1) / 2)),
      m_mapToPair(m_distances.size()),
      m_activeToOriginal(cmps.size()),
      m_originalToActive(cmps.size()),
      m_numberComponents(cmps.size()),
      m_numberActive(cmps.size()) {
  std::iota(m_activeToOriginal.begin(), m_activeToOriginal.end(), 0);
  std::iota(m_originalToActive.begin(), m_originalToActive.end(), 0);

  for (std::size_t i = 1; i < m_numberComponents; ++i) {
    const std::size_t indexConst = (i - 1) * i / 2;
    for (std::size_t j = 0; j < i; ++j) {
      m_mapToPair.at(indexConst + j) = {i, j};
      m_distances[indexConst + j] =
          computeSymmetricKlDivergence(cmps[i], cmps[j]);
    }
  }
}

void SymmetricKLDistanceMatrix::recomputeAssociatedDistances(
    std::size_t originalIdx, std::span<const GsfComponent> cmps) {
  assert(cmps.size() == m_numberComponents && "size mismatch");

  const std::size_t activeIdx = m_originalToActive[originalIdx];

  setAssociated(activeIdx, [&cmps, this](std::size_t row, std::size_t col) {
    return computeSymmetricKlDivergence(cmps[m_activeToOriginal[row]],
                                        cmps[m_activeToOriginal[col]]);
  });
}

void SymmetricKLDistanceMatrix::maskAssociatedDistances(
    std::size_t originalIdx) {
  const std::size_t activeRemoved = m_originalToActive[originalIdx];
  const std::size_t last = m_numberActive - 1;

  if (activeRemoved != last) {
    // Rows in distance matrix
    for (std::size_t i = 0; i < activeRemoved; ++i) {
      std::swap(m_distances[triangularIndex(activeRemoved, i)],
                m_distances[triangularIndex(last, i)]);
    }
    // Columns in distance matrix
    for (std::size_t i = activeRemoved + 1; i < last; ++i) {
      std::swap(m_distances[triangularIndex(i, activeRemoved)],
                m_distances[triangularIndex(last, i)]);
    }

    // Move the last active component into the freed slot
    const std::size_t lastOriginal = m_activeToOriginal[last];
    m_activeToOriginal[activeRemoved] = lastOriginal;
    m_originalToActive[lastOriginal] = activeRemoved;
  }

  m_numberActive = last;
}

std::pair<std::size_t, std::size_t> SymmetricKLDistanceMatrix::minDistancePair()
    const {
  const std::size_t nActivePairs = m_numberActive * (m_numberActive - 1) / 2;

  // Single-pass, block-tracking argmin: scan in fixed-size blocks with a
  // plain sequential min reduction (auto-vectorizable, since min is
  // associative/commutative), tracking which block held the running
  // minimum. Because the minimum only updates on a strict '<', that block
  // is guaranteed to contain the first occurrence, so recovering the index
  // only needs a small rescan of one block instead of the whole array.
  // BLOCK=16 was empirically best for the N=48-72 active components (up to
  // ~2556 pairs) typical here.
  constexpr std::size_t BLOCK = 16;

  const Array &data = m_distances;
  double currentMin = std::numeric_limits<double>::max();
  std::size_t currentBlockStart = 0;
  std::size_t currentBlockLen = 0;

  std::size_t i = 0;
  for (; i + BLOCK <= nActivePairs; i += BLOCK) {
    double blockMin = data[i];
    for (std::size_t k = 1; k < BLOCK; ++k) {
      blockMin = std::min(blockMin, data[i + k]);
    }
    if (blockMin < currentMin) {
      currentMin = blockMin;
      currentBlockStart = i;
      currentBlockLen = BLOCK;
    }
  }
  if (i < nActivePairs) {
    double blockMin = data[i];
    for (std::size_t k = i + 1; k < nActivePairs; ++k) {
      blockMin = std::min(blockMin, data[k]);
    }
    if (blockMin < currentMin) {
      currentMin = blockMin;
      currentBlockStart = i;
      currentBlockLen = nActivePairs - i;
    }
  }

  // Find the first occurrence of the minimum within the winning block; tie
  // order doesn't matter since merging exactly-tied components is
  // equivalent either way.
  std::size_t idx = currentBlockStart;
  for (std::size_t j = currentBlockStart;
       j < currentBlockStart + currentBlockLen; ++j) {
    if (data[j] == currentMin) {
      idx = j;
      break;
    }
  }

  const auto [row, col] = m_mapToPair[idx];
  return {m_activeToOriginal[row], m_activeToOriginal[col]};
}

std::ostream &SymmetricKLDistanceMatrix::toStream(std::ostream &os) const {
  const int width = 8;
  const int prec = 2;

  os << "\n";
  os << std::format("{:>{}}", "", width) << " | ";
  for (std::size_t j = 0ul; j < m_numberActive - 1; ++j) {
    os << std::format("{:>{}}", m_activeToOriginal[j], width) << "  ";
  }
  os << "\n";
  os << std::string((width + 3) + (width + 2) * (m_numberActive - 1), '-');
  os << "\n";

  for (std::size_t i = 1ul; i < m_numberActive; ++i) {
    os << std::format("{:>{}}", m_activeToOriginal[i], width) << " | ";
    for (std::size_t j = 0ul; j < i; ++j) {
      os << std::format("{:>{}.{}f}", m_distances[triangularIndex(i, j)], width,
                        prec)
         << "  ";
    }
    os << "\n";
  }
  return os;
}

}  // namespace detail::Gsf

}  // namespace Acts

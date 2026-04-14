// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/CombinatoricIndices.hpp"

#include <format>

namespace Acts {

template <std::size_t K>
CombinatoricIndices<K>::CombinatoricIndices(const std::size_t N) : m_N{N} {
  if (N < K) {
    throw std::invalid_argument(
        std::format("CombinatoricIndices() - The set size {:} needs at least "
                    "to exceed  the number of elements to draw {:}",
                    N, K));
  }
  static_assert(K >= 1, "The number of elements must not be zero");

  /// Use the identity (N, M) = (N, N-M)
  ///  (N, K) = sum_{i=K)^{N-1} (I, K)
  /// Use the identity (N, M) = (N, N-M)
  ///  (N, K) = sum_{i=K)^{N-1} (I, K)
  m_borders.reserve(N - K);
  std::size_t setSize{0ul};
  for (std::size_t i = N - 1; i >= K - 1ul; --i) {
    setSize += binomial(i, K - 1);
    m_borders.push_back(setSize);
    if (i == 0) {
      break;
    }
  }
}

template <std::size_t K>
std::size_t CombinatoricIndices<K>::size() const {
  return m_borders.back();
}
template <std::size_t K>
std::size_t CombinatoricIndices<K>::setSize() const {
  return m_N;
}

template <std::size_t K>
std::array<std::size_t, K> CombinatoricIndices<K>::draw(
    const std::size_t combination) const {
  if (combination >= size()) {
    throw std::invalid_argument(
        std::format("CombinatoricIndices::draw({:}) - You cannot draw more "
                    "than {:} combinations",
                    combination, size()));
  }
  std::array<std::size_t, K> indices{};
  for (std::size_t s = 0; s < indices.size(); ++s) {
    indices[s] = drawIndex(combination, s);
  }
  return indices;
}

template <std::size_t K>

std::size_t CombinatoricIndices<K>::drawIndex(const std::size_t combination,
                                              const std::size_t slot) const {
  // There is only one combination if N == K
  if (m_N == K) {
    return slot;
  }
  if (slot >= K) {
    throw std::invalid_argument(std::format(
        "CombinatoricIndices - The slot {:} must not exceed {:}.", slot, K));
  }
  const auto slotItr = std::ranges::upper_bound(m_borders, combination);
  if (slot == 0ul) {
    return std::distance(m_borders.begin(), slotItr);
  }
  if constexpr (s_hasSubDraw) {
    const std::size_t combInSubDraw =
        m_borders.front() - ((*slotItr) - combination);
    return m_childDrawer.drawIndex(combInSubDraw, slot - 1ul) + 1ul;
  } else {
    throw std::invalid_argument(
        std::format("CombinatoricIndices ({:}, {:}). The combination {:} "
                    "exceeds the allowed range {:}",
                    m_N, K, combination, size()));
    return 0ul;
  }
}

}  // namespace Acts

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/detail/CombinatoricIndices.hpp"

#include <format>

namespace Acts::detail {

template <std::size_t K>
CombinatoricIndices<K>::CombinatoricIndices(const std::size_t N) : m_N{N} {
  if (N < K) {
    throw std::invalid_argument(
        std::format("CombinatoricIndices() - The set size {:} needs at least "
                    "to exceed the number of elements to draw {:}",
                    N, K));
  }
  static_assert(K >= 1, "The number of elements must not be zero");

  /// Use the identity (N, M) = (N, N-M)
  ///  (N, K) = sum_{i=K)^{N-1} (I, K)
  /// Use the identity (N, M) = (N, N-M)
  ///  (N, K) = sum_{i=K)^{N-1} (I, K)
  m_borders.reserve(N - K);
  std::size_t setSize{0ul};
  /// Calculate the number of combinations in which the i-th element
  /// in the sequence is used
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
    const std::size_t subCombination =
        m_borders.front() - ((*slotItr) - combination);
    return m_childDrawer.drawIndex(subCombination, slot - 1ul) + 1ul;
  } else {
    throw std::invalid_argument(
        std::format("CombinatoricIndices ({:}, {:}). The combination {:} "
                    "exceeds the allowed range {:}",
                    m_N, K, combination, size()));
    return 0ul;
  }
}

template <std::size_t K>
CombinatoricIndices<K>::iterator::iterator(const CombinatoricIndices* parent,
                                           const std::size_t _itr)
    : m_parent{parent}, m_itr{_itr} {
  updateArray();
}

template <std::size_t K>
CombinatoricIndices<K>::iterator&
CombinatoricIndices<K>::iterator::operator++() {
  ++m_itr;
  updateArray();
  return *this;
}

template <std::size_t K>
bool CombinatoricIndices<K>::iterator::operator==(const iterator& other) const {
  return m_parent == other.m_parent && m_itr == other.m_itr;
}

template <std::size_t K>
bool CombinatoricIndices<K>::iterator::operator!=(const iterator& other) const {
  return m_parent != other.m_parent || m_itr != other.m_itr;
}

template <std::size_t K>
const CombinatoricIndices<K>::IndexArray&
CombinatoricIndices<K>::iterator::operator*() const {
  return m_array;
}

template <std::size_t K>
void CombinatoricIndices<K>::iterator::updateArray() {
  if (m_parent == nullptr || m_itr >= m_parent->size()) {
    m_array =
        filledArray<std::size_t, K>(std::numeric_limits<std::size_t>::max());
  } else {
    for (std::size_t s = 0; s < m_array.size(); ++s) {
      m_array[s] = m_parent->drawIndex(m_itr, s);
    }
  }
}
template <std::size_t K>
CombinatoricIndices<K>::iterator CombinatoricIndices<K>::iterator::operator+(
    const std::size_t idx) const {
  return iterator{m_parent, m_itr + idx};
}

}  // namespace Acts::detail

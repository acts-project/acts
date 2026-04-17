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
CombinatoricIndices<K>::CombinatoricIndices(const std::size_t N)
  requires(K > 0)
    : m_N{N} {
  if (N < K) {
    throw std::invalid_argument(
        std::format("CombinatoricIndices() - The set size {:} needs at least "
                    "to exceed the number of elements to draw {:}",
                    N, K));
  }

  for (std::size_t slot = 0; slot < m_elementOccurance.size(); ++slot) {
    std::size_t maxCombNumb = 0;
    // Calculate in how many combinations the integer can occur
    const std::size_t nPrime = N - 1ul - slot;
    const std::size_t kPrime = K - 1ul - slot;

    for (std::size_t leftCmb = nPrime; leftCmb >= kPrime; --leftCmb) {
      maxCombNumb += binomial(leftCmb, kPrime);
      m_elementOccurance[slot].push_back(maxCombNumb);
      if (leftCmb == 0ul) {
        break;
      }
    }
  }
}

template <std::size_t K>
std::size_t CombinatoricIndices<K>::size() const {
  return m_elementOccurance[0].back();
}

template <std::size_t K>
std::size_t CombinatoricIndices<K>::intervalSize() const {
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
  std::size_t cmbPrime = combination;
  if (slot >= K) {
    throw std::runtime_error(
        std::format("CombinatoricIndices: The slot must not exceed {:} the "
                    "number of drawn elements {:}",
                    slot, K));
  }
  for (std::size_t s = 0; s <= slot; ++s) {
    const auto slotItr =
        std::ranges::upper_bound(m_elementOccurance[s], cmbPrime);
    if (s == slot) {
      return std::distance(m_elementOccurance[s].begin(), slotItr) + slot;
    }
    cmbPrime -= ((*slotItr) - m_elementOccurance[s].front());
  }
  return 0ul;
}

template <std::size_t K>
CombinatoricIndices<K>::iterator::iterator(const CombinatoricIndices* parent,
                                           const std::size_t combination)
    : m_parent{parent}, m_itr{combination} {
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

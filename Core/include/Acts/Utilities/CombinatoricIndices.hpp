// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/ArrayHelpers.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <vector>

namespace Acts {
/// Utilitty class to draw
template <std::size_t K>
class CombinatoricIndices {
 public:
  /// Declare other combinatoric indices classes as friends
  template <std::size_t L>
  friend class CombinatoricIndices;

  using IndexArray = std::array<std::size_t, K>;

  /// Constructor of the combinatoric indices
  /// @param The size of the set from which indices are drawn
  /// @note An exception is thrown if the size is less than the
  /// number of indices to draw
  explicit CombinatoricIndices(const std::size_t N);

  /// The number of possibilities to draw K different inidces from the set
  /// @returns The number of combinations
  std::size_t size() const;
  /// The size of the set from which the indices are drawn
  /// @returns The parameter N with which this class was constructed
  std::size_t setSize() const;

  /// Draws a new combination of indices and stores them in the array
  /// @param combination: The number of the combination in the sequence
  /// @returns An array where each is a unique number from [0 -N)
  IndexArray draw(const std::size_t combination) const;

  class iterator {
   public:
    iterator() = default;
    iterator(const CombinatoricIndices* parent, const std::size_t _itr)
        : m_parent{parent}, m_itr{_itr} {
      updateArray();
    }

    iterator& operator++() {
      ++m_itr;
      updateArray();
      return *this;
    }
    bool operator==(const iterator& other) const {
      return m_parent == other.m_parent && m_itr == other.m_itr;
    }
    bool operator!=(const iterator& other) const {
      return m_parent != other.m_parent || m_itr != other.m_itr;
    }
    const IndexArray& operator*() const { return m_array; }

   private:
    void updateArray() {
      if (m_parent == nullptr || m_itr >= m_parent->size()) {
        m_array = filledArray<std::size_t, K>(
            std::numeric_limits<std::size_t>::max());
      } else {
        for (std::size_t s = 0; s < m_array.size(); ++s) {
          m_array[s] = m_parent->drawIndex(m_itr, s);
        }
      }
    }
    const CombinatoricIndices* m_parent{nullptr};
    std::size_t m_itr{0ul};
    IndexArray m_array{
        filledArray<std::size_t, K>(std::numeric_limits<std::size_t>::max())};
  };

  iterator begin() const { return iterator{this, 0ul}; }

  iterator end() const { return iterator{this, size()}; }

 private:
  std::size_t drawIndex(const std::size_t combination,
                        const std::size_t slot) const;

  static constexpr bool s_hasSubDraw = K > 1ul;
  using SubDrawer_t =
      std::conditional_t<s_hasSubDraw, CombinatoricIndices<K - 1>, std::size_t>;

  std::size_t m_N{};
  SubDrawer_t m_childDrawer{m_N - 1ul};

  std::vector<std::size_t> m_borders{};
};

}  // namespace Acts

#include "Acts/Utilities/CombinatoricIndices.ipp"

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

namespace Acts::detail {
/// Utility class to loop over all possible ways to draw K unique elements out
/// of a continuous sequence of N elements. The list of combinations starts with
/// the lowest possible set of indices, e.g. for a sequence of 4
///             [0, 1, 2, 3]
/// The next members of the series increment the (K-1) th index until
///             [0, 1, 2, N-1].
/// Then the third index is incremented and the last one set to the next value
///             [0, 1, 3, 4].
/// The sequence continues until all elements are at their maximum value.
///
/// The user can either loop manually to draw the list of elements
/// or use a range based for loop
template <std::size_t K>
class CombinatoricIndices {
 public:
  /// Declare the Return type of the Index generator to be an array of size K
  using IndexArray = std::array<std::size_t, K>;

  /// Delete default constructor
  CombinatoricIndices() = delete;
  /// Constructor of the combinatoric indices
  /// @param n The size of the set from which indices are drawn
  /// @note An exception is thrown if the size is less than the
  /// number of indices to draw
  explicit CombinatoricIndices(const std::size_t N)
    requires(K > 0);

  /// The number of possibilities to draw K different indices from the set
  /// @returns The number of combinations
  std::size_t size() const;
  /// The size of the set from which the indices are drawn
  /// @returns The parameter N with which this class was constructed
  std::size_t intervalSize() const;

  /// Draws a new combination of indices and stores them in the array
  /// @param combination The number of the combination in the sequence
  /// @returns An array where each is a unique number from [0; N)
  IndexArray draw(const std::size_t combination) const;

  /// Iterator class to be used in the ranged for loop
  class iterator {
   public:
    /// Empty default constructor
    iterator() = default;
    /// Constructor with a CombinatoricIndex parent and the iterator position
    /// @param parent Pointer to the parent used to draw the combinatoric indices
    /// @param combination Position of the iterator in the sequence of combinations
    iterator(const CombinatoricIndices* parent, const std::size_t combination);

    /// Increment the internal iterator count by one unit
    /// @returns Mutable reference to the iterator instance
    iterator& operator++();
    /// Return an iterator shifted by x indices in the sequence
    /// @param idx Number of indexes by which the new operator is shifted from this one
    /// @returns A new iterator
    iterator operator+(const std::size_t idx) const;

    /// Comparison operator to check whether two iterators are equal,
    /// i.e. they share the same parent and have the same internal iterator
    /// index
    /// @param other Const reference to the other iterator to check
    /// @returns The equality assessment of the two iterators
    bool operator==(const iterator& other) const;

    /// Comparison operator to check whether two iterators are unequal
    /// @param other Const reference to the other iterator to check
    /// @returns The inequality assessment of the two iterators
    bool operator!=(const iterator& other) const;

    /// Dereference operator to the underlying memory
    /// @returns The array containing the unique indices of the combinatoric generator
    const IndexArray& operator*() const;

   private:
    /// Update the internal array to the current element in the sequence
    void updateArray();
    /// Pointer to the parent from which the combinatoric indices are drawn
    const CombinatoricIndices* m_parent{nullptr};
    /// Iterator index in the sequence of iterations
    std::size_t m_itr{0ul};
    /// Storage for the drawn indices
    IndexArray m_array{
        filledArray<std::size_t, K>(std::numeric_limits<std::size_t>::max())};
  };

  /// @returns the begin iterator of the sequence of combinatoric indices
  constexpr iterator begin() const { return iterator{this, 0ul}; }
  /// @returns the end iterator of the sequence of combinatoric indices
  constexpr iterator end() const { return iterator{this, size()}; }

 private:
  /// Draws the unique index
  /// @param combination: The number of the combination in the list
  ///         of sequences
  /// @param slot The positional index of the index inside the returned array
  /// @returns A unique index for the n-th combination
  std::size_t drawIndex(const std::size_t combination,
                        const std::size_t slot) const;

  /// Size of the set from which the combinations can be drawn needs to be >= K
  std::size_t m_N{};
  /// Cache to state until which member in the sequence of drawn combinations
  /// the lowest possible integer can occur in a slot. E.g. for the case N = 5,
  /// k = 3 The 0 at the first position is used until combination #6, the 1
  /// until combination 6 + 3 = 9 and the 2 until combination  6 + 3 + 1 = 10
  std::array<std::vector<std::size_t>, K> m_elementOccurance{};
};

}  // namespace Acts::detail

#include "Acts/Utilities/detail/CombinatoricIndices.ipp"

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>

namespace Acts {

/// This object can be iterated to produce up to two sequences of integer
/// indices, corresponding to the half-open integer ranges [begin1, end1[ and
/// [begin2, end2[.
///
/// The goal is to emulate the effect of enumerating a range of neighbor
/// indices on an axis (which may go out of bounds and wrap around since we
/// have AxisBoundaryType::Closed), inserting them into an std::vector, and
/// discarding duplicates, without paying the price of duplicate removal
/// and dynamic memory allocation in hot magnetic field interpolation code.
///
/// Iterable indices for neighborhood lookups with optional wrap-around.
class NeighborHoodIndices {
 public:
  NeighborHoodIndices() = default;

  /// Constructor for continuous range
  /// @param begin Start index
  /// @param end End index (exclusive)
  NeighborHoodIndices(std::size_t begin, std::size_t end)
      : m_begin1(begin), m_end1(end), m_begin2(end), m_end2(end) {}

  /// Constructor for wrapped range (two segments)
  /// @param begin1 Start of first segment
  /// @param end1 End of first segment (exclusive)
  /// @param begin2 Start of second segment
  /// @param end2 End of second segment (exclusive)
  NeighborHoodIndices(std::size_t begin1, std::size_t end1, std::size_t begin2,
                      std::size_t end2)
      : m_begin1(begin1), m_end1(end1), m_begin2(begin2), m_end2(end2) {}

  /// Iterator over the neighborhood index sequence.
  class iterator {
   public:
    iterator() = default;

    /// Constructor for end iterator
    /// @param current End position
    explicit iterator(std::size_t current)
        : m_current(current), m_wrapped(true) {}

    /// Constructor for begin iterator
    /// @param begin1 Start of first segment
    /// @param end1 End of first segment
    /// @param begin2 Start of second segment
    iterator(std::size_t begin1, std::size_t end1, std::size_t begin2)
        : m_current(begin1),
          m_end1(end1),
          m_begin2(begin2),
          m_wrapped(begin1 == begin2) {}

    /// Dereference operator
    /// @return Current index
    std::size_t operator*() const { return m_current; }

    /// Pre-increment operator
    /// @return Reference to this iterator
    iterator& operator++() {
      ++m_current;
      if (m_current == m_end1) {
        m_current = m_begin2;
        m_wrapped = true;
      }
      return *this;
    }

    /// Equality comparison operator
    /// @param it Other iterator
    /// @return True if iterators are equal
    bool operator==(const iterator& it) const {
      return (m_current == it.m_current) && (m_wrapped == it.m_wrapped);
    }

   private:
    std::size_t m_current = 0, m_end1 = 0, m_begin2 = 0;
    bool m_wrapped = false;
  };

  /// Get begin iterator
  /// @return Iterator to first index
  iterator begin() const { return iterator(m_begin1, m_end1, m_begin2); }

  /// Get end iterator
  /// @return Iterator past last index
  iterator end() const { return iterator(m_end2); }

  /// Get total number of indices in the sequence
  /// @return Number of indices
  std::size_t size() const { return (m_end1 - m_begin1) + (m_end2 - m_begin2); }

  /// Collect all indices into a vector
  /// @return Vector containing all indices
  std::vector<std::size_t> collect() const {
    std::vector<std::size_t> result;
    result.reserve(this->size());
    for (std::size_t idx : *this) {
      result.push_back(idx);
    }
    return result;
  }

 private:
  std::size_t m_begin1 = 0;
  std::size_t m_end1 = 0;
  std::size_t m_begin2 = 0;
  std::size_t m_end2 = 0;
};

}  // namespace Acts

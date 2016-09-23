// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_SEEDING_RANGE_HPP
#define ACTS_SEEDING_RANGE_HPP

#include <algorithm>

#include "ACTS/Utilities/detail/periodic.hpp"

namespace Acts {
namespace detail {

  /// A range of elements inside a container that wraps around the end.
  ///
  /// @note Container must be indexable, i.e. provide `.size()` and
  ///       `.operator[](...)`.
  template <typename Container, typename Index = typename Container::size_type>
  struct CyclicRange
  {
    struct const_iterator
    {
      using difference_type   = typename Container::difference_type;
      using value_type        = typename Container::value_type;
      using pointer           = typename Container::pointer;
      using reference         = typename Container::reference;
      using iterator_category = std::forward_iterator_tag;

      const Container& values;
      Index            index;

      const value_type& operator*() const
      {
        return values[index % values.size()];
      }
      bool
      operator!=(const const_iterator& other) const
      {
        return (index != other.index);
      }
      const_iterator& operator++()
      {
        ++index;
        return *this;
      }
    };

    const Container& values;
    Index            start, length;

    const_iterator
    begin() const
    {
      return {values, start};
    }
    const_iterator
    end() const
    {
      return {values, start + length};
    }
    bool
    empty() const
    {
      return (length == 0);
    }
  };

  /// Select subset of ordered points within phi range.
  ///
  /// The container must fulfill the preconditions for the `CyclicRange`.
  /// In addition, its iterators must be RandomAccessIterators, its values must
  /// provide a `.phi()` accessor that returns values in the [-phi, phi) range
  /// and the container must be ordered with respect to this value.
  template <typename Container, typename Index = typename Container::size_type>
  CyclicRange<Container, Index>
  makeRangePhi(const Container& values, double phi0, double phi1)
  {
    using Value = typename Container::value_type;

    phi0 = Acts::detail::radian_sym(phi0);
    phi1 = Acts::detail::radian_sym(phi1);

    // assumes a sorted container for fast upper/lower binary search
    auto cmpValuePhi = [](const Value& a, double b) { return (a.phi() < b); };
    auto cmpPhiValue = [](double a, const Value& b) { return (a < b.phi()); };
    auto it0
        = std::lower_bound(values.begin(), values.end(), phi0, cmpValuePhi);
    // all elements are smaller then lower limit -> empty range
    if (it0 == values.end()) {
      return {values, 0, 0};
    }
    auto it1
        = std::upper_bound(values.begin(), values.end(), phi1, cmpPhiValue);
    auto begin = std::distance(values.begin(), it0);
    auto diff  = std::distance(it0, it1);
    if (diff < 0) {
      // range wraps around the end
      return {values, Index(begin), Index(values.size() + diff)};
    }
    return {values, Index(begin), Index(diff)};
  }

}  // namespace detail
}  // namespace Acts

#endif  // ACTS_SEEDING_RANGE_HPP

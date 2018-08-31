// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <algorithm>

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

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

  template <typename Container, typename Index>
  std::ostream&
  operator<<(std::ostream& os, const CyclicRange<Container, Index>& range)
  {
    if (range.empty()) {
      os << '{' << range.start << ',' << range.start << '}';
    } else {
      auto first = range.start;
      auto last  = range.start + range.length - 1;
      os << '{' << first << "=(" << range.values[first] << "), " << last << "=("
         << range.values[last % range.values.size()] << ")}";
    }
    return os;
  }

  /// Select subset of ordered points within phi range.
  ///
  /// @note Assumes a sorted container for fast upper/lower binary search.
  ///
  /// The container must fulfill the preconditions for the `CyclicRange`.
  /// In addition, its iterators must be RandomAccessIterators, its values must
  /// provide a `.phi()` accessor that returns values in the [-phi, phi) range
  /// and the container must be ordered with respect to this value.
  template <typename Container, typename Index = typename Container::size_type>
  CyclicRange<Container, Index>
  makeRangePhi(const Container& values, double phi0, double phi1)
  {
    using VectorHelpers::phi;
    phi0 = Acts::detail::radian_sym(phi0);
    phi1 = Acts::detail::radian_sym(phi1);

    auto compValPhi = [](const auto& a, double b) { return (phi(a) < b); };
    auto compPhiVal = [](double a, const auto& b) { return (a < phi(b)); };

    // linear boundaries that do not allow wrap-around
    if (phi0 <= phi1) {
      auto it0
          = std::lower_bound(values.begin(), values.end(), phi0, compValPhi);
      auto it1 = std::upper_bound(it0, values.end(), phi1, compPhiVal);
      // special value for empty range
      if (it0 == it1) {
        return {values, 0, 0};
      }
      return {values,
              Index(std::distance(values.begin(), it0)),
              Index(std::distance(it0, it1))};
    }
    // inverted boundaries w/ possible wrap-around
    auto it1 = std::upper_bound(values.begin(), values.end(), phi1, compPhiVal);
    auto it0 = std::lower_bound(it1, values.end(), phi0, compValPhi);
    // wrap-around can be reduced to linear part
    if (it0 == values.end()) {
      return {values, 0, Index(std::distance(values.begin(), it1))};
    }
    // full wrap-around w/ elements before and after end
    return {values,
            Index(std::distance(values.begin(), it0)),
            Index(values.size() - std::distance(it1, it0))};
  }

}  // namespace detail
}  // namespace Acts
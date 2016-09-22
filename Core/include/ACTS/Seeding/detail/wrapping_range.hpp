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
namespace Seeding {
  namespace detail {

    /// An iterator that wraps at the end of the container.
    ///
    /// The container must be indexable, i.e. must provide a `size()` method and
    /// an `operator[]` to access elements.
    template <typename Container,
              typename Index = typename Container::size_type>
    struct WrappingConstIterator
    {
      const Container& values;
      Index            index;

      auto operator*() const { return values[index % values.size()]; }
      auto operator-> () const { return &values[index % values.size()]; }
      bool
      operator!=(const WrappingConstIterator& other) const
      {
        return index != other.index;
      }
      bool
      operator==(const WrappingConstIterator& other) const
      {
        return index == other.index;
      }
      WrappingConstIterator& operator++()
      {
        ++index;
        return *this;
      }
    };

    /// A range of elements inside a container that can wrap around the end.
    ///
    /// The container must be indexable, i.e. must provide a `size()` method
    /// and an `operator[]` to access elements.
    template <typename Container,
              typename Index = typename Container::size_type>
    struct WrappingRange
    {
      using Iterator = WrappingConstIterator<Container, Index>;

      const Container& values;
      Index            start, length;

      Iterator
      begin() const
      {
        return {values, start};
      }
      Iterator
      end() const
      {
        return {values, start + length};
      }
      bool
      empty() const
      {
        return (length == 0);
      }
      Index
      size() const
      {
        return length;
      }
    };

    /// Select subset of ordered points within phi range.
    ///
    /// The container must fulfill the preconditions for the `WrappingRange`.
    /// In addition, its values must provide a `.phi()` accessor that returns
    /// values in the [-phi, phi) range and the container must be ordered with
    /// respect to this value.
    template <typename Container,
              typename Index = typename Container::size_type,
              typename Value = typename Container::value_type>
    WrappingRange<Container, Index>
    makeRangePhi(const Container& values, double phi0, double phi1)
    {
      phi0 = Acts::detail::radian_sym(phi0);
      phi1 = Acts::detail::radian_sym(phi1);

      // std algorithms assumes a sorted container for fast upper/lower
      // bound binary search.
      auto cmpValuePhi = [](const Value& a, double b) { return (a.phi() < b); };
      auto cmpPhiValue = [](double a, const Value& b) { return (a < b.phi()); };
      auto it0
          = std::lower_bound(values.begin(), values.end(), phi0, cmpValuePhi);
      auto it1
          = std::upper_bound(values.begin(), values.end(), phi1, cmpPhiValue);

      // all elements are smaller then lower limit -> empty range
      if (it0 == values.end()) return {values, 0, 0};

      size_t start = std::distance(values.begin(), it0);

      // no elements are larger than the upper limit -> remainder of values
      if (it1 == values.end()) return {values, start, values.size() - start};

      size_t upper = std::distance(values.begin(), it1);
      // no cyclic wrapping, range is fully linear
      // NOTE: upper bound is first element not in limits
      if (start <= upper) return {values, start, upper - start};

      // selected range wraps at the end
      return {values, start, values.size() - start + upper};
    }

  }  // namespace detail
}  // namespace Seeding
}  // namespace Acts

#endif  // ACTS_SEEDING_RANGE_HPP

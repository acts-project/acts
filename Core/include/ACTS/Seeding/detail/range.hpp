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

    /**
     * @brief An iterator that automatically wraps at the end of the container.
     *
     * The container must be indexable, i.e. must provide a `size()` method and
     * an `operator[]` to access elements.
     */
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

    /**
     * @brief An subset of ordered points within some phi range.
     *
     * The container must be indexable, i.e. must provide a `size()` method and
     * an `operator[]` to access elements. The stored value must provide a
     * `phi()` method that returns the azimuthal angle in the [-pi,pi) range.
     */
    template <typename Container,
              typename Index = typename Container::size_type,
              typename Value = typename Container::value_type>
    struct PhiRange
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

      /**
       * @brief Select elements in the half-open [phi0,phi1) interval.
       *
       * Periodicity in phi is automatically considered.
       */
      PhiRange(const Container& values_, double phi0, double phi1)
        : values(values_)
      {
        phi0 = Acts::detail::radian_sym(phi0);
        phi1 = Acts::detail::radian_sym(phi1);

        // std algorithms assumes a sorted container for fast upper/lower
        // bound binary search.
        auto cmpValuePhi
            = [](const Value& a, double b) { return (a.phi() < b); };
        auto cmpPhiValue
            = [](double a, const Value& b) { return (a < b.phi()); };
        auto it0 = std::lower_bound(
            values_.begin(), values_.end(), phi0, cmpValuePhi);
        auto it1 = std::upper_bound(
            values_.begin(), values_.end(), phi1, cmpPhiValue);

        if (it0 == values_.end()) {
          // all elements are smaller then lower limit -> empty range
          start  = 0;
          length = 0;
        } else {
          start = std::distance(values_.begin(), it0);
          // NOTE: upper bound is first element not in limits
          // NOTE: upper bound can also be values_.end()
          auto upper = std::distance(values_.begin(), it1);
          if (start <= upper) {
            // no cyclic wrapping, range is fully linear
            length = upper - start;
          } else {
            // partial range at the end of the container
            length = values_.size() - start;
            // partial range wrapped back to the beginning
            length += upper;
          }
        }
      }
    };

    template <typename Container>
    PhiRange<Container>
    makeRangePhi(const Container& values, double phi0, double phi1)
    {
      return {values, phi0, phi1};
    }

  }  // namespace detail
}  // namespace Seeding
}  // namespace Acts

#endif  // ACTS_SEEDING_RANGE_HPP

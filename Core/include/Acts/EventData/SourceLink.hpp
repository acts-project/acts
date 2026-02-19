// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Any.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <cassert>
#include <concepts>
#include <type_traits>
#include <utility>

#if !defined(ACTS_SOURCELINK_SBO_SIZE)
// Do not set this in code, use CMake!
#define ACTS_SOURCELINK_SBO_SIZE 16
#endif

namespace Acts {

/// @ingroup eventdata_measurement
/// Type-erased source link wrapper.
class SourceLink final {
  using any_type = AnyBase<ACTS_SOURCELINK_SBO_SIZE>;

 public:
  /// Copy constructor
  /// @param other Source link to copy from
  SourceLink(const SourceLink& other) = default;
  /// Move constructor
  /// @param other Source link to move from
  SourceLink(SourceLink&& other) = default;
  /// Copy assignment operator
  /// @param other Source link to copy from
  /// @return Reference to this source link
  SourceLink& operator=(const SourceLink& other) = default;
  /// Move assignment operator
  /// @param other Source link to move from
  /// @return Reference to this source link
  SourceLink& operator=(SourceLink&& other) = default;

  /// Constructor from concrete source link
  /// @tparam T The source link type
  /// @param upstream The upstream source link to store
  template <typename T>
  explicit SourceLink(T&& upstream)
    requires(!std::same_as<std::decay_t<T>, SourceLink>)
      : m_upstream(std::forward<T>(upstream)) {
    static_assert(!std::is_same_v<std::decay_t<T>, SourceLink>,
                  "Cannot wrap SourceLink in SourceLink");
  }

  /// Concrete source link class getter
  /// @tparam T The source link type to retrieve
  /// @return Reference to the stored source link
  template <typename T>
  T& get() {
    return m_upstream.as<T>();
  }

  /// Concrete source link class getter, const version
  /// @tparam T The source link type to retrieve
  /// @return Const reference to the stored source link
  template <typename T>
  const T& get() const {
    return m_upstream.as<T>();
  }

 private:
  any_type m_upstream{};
};

/// Iterator adapter returning SourceLink wrappers.
template <typename T>
struct SourceLinkAdapterIterator {
  /// The base iterator type
  using BaseIterator = T;

  /// Iterator category
  using iterator_category = typename BaseIterator::iterator_category;
  /// Value type
  using value_type = typename BaseIterator::value_type;
  /// Difference type
  using difference_type = typename BaseIterator::difference_type;
  /// Pointer type
  using pointer = typename BaseIterator::pointer;
  /// Reference type
  using reference = typename BaseIterator::reference;

  /// Constructor from base iterator
  /// @param iterator Base iterator to wrap
  explicit SourceLinkAdapterIterator(const T& iterator)
      : m_iterator{iterator} {}

  /// Pre-increment operator
  /// @return Reference to this iterator
  SourceLinkAdapterIterator& operator++() {
    ++m_iterator;
    return *this;
  }

  /// Equality comparison operator
  /// @param other Iterator to compare with
  /// @return True if iterators are equal
  bool operator==(const SourceLinkAdapterIterator& other) const {
    return m_iterator == other.m_iterator;
  }

  /// Dereference operator
  /// @return Wrapped source link
  Acts::SourceLink operator*() const { return Acts::SourceLink{*m_iterator}; }

  /// Difference operator
  /// @param other Iterator to compute difference with
  /// @return Distance between iterators
  auto operator-(const SourceLinkAdapterIterator& other) const {
    return m_iterator - other.m_iterator;
  }

  /// Underlying iterator
  BaseIterator m_iterator;
};

/// Delegate to unpack the surface associated with a source link
using SourceLinkSurfaceAccessor = Delegate<const Surface*(const SourceLink&)>;

}  // namespace Acts

// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Any.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <cassert>
#include <iostream>
#include <type_traits>
#include <utility>

#if !defined(ACTS_SOURCELINK_SBO_SIZE)
// Do not set this in code, use CMake!
#define ACTS_SOURCELINK_SBO_SIZE 16
#endif

namespace Acts {

class SourceLink final {
  using any_type = AnyBase<ACTS_SOURCELINK_SBO_SIZE>;

 public:
  SourceLink(const SourceLink& other) = default;
  SourceLink(SourceLink&& other) = default;
  SourceLink& operator=(const SourceLink& other) = default;
  SourceLink& operator=(SourceLink&& other) = default;

  /// Constructor from concrete sourcelink
  /// @tparam T The source link type
  /// @param upstream The upstream source link to store
  template <typename T, typename = std::enable_if_t<
                            !std::is_same_v<std::decay_t<T>, SourceLink>>>
  explicit SourceLink(T&& upstream) {
    static_assert(!std::is_same_v<std::decay_t<T>, SourceLink>,
                  "Cannot wrap SourceLink in SourceLink");

    if constexpr (std::is_same_v<T, std::decay_t<T>>) {
      m_upstream = any_type{std::move(upstream)};
    } else {
      m_upstream = any_type{static_cast<std::decay_t<T>>(upstream)};
    }
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

template <typename T>
struct SourceLinkAdapterIterator {
  using BaseIterator = T;

  using iterator_category = typename BaseIterator::iterator_category;
  using value_type = typename BaseIterator::value_type;
  using difference_type = typename BaseIterator::difference_type;
  using pointer = typename BaseIterator::pointer;
  using reference = typename BaseIterator::reference;

  explicit SourceLinkAdapterIterator(T iterator) : m_iterator{iterator} {}

  SourceLinkAdapterIterator& operator++() {
    ++m_iterator;
    return *this;
  }

  bool operator==(const SourceLinkAdapterIterator& other) const {
    return m_iterator == other.m_iterator;
  }

  bool operator!=(const SourceLinkAdapterIterator& other) const {
    return !(*this == other);
  }

  Acts::SourceLink operator*() const { return Acts::SourceLink{*m_iterator}; }

  auto operator-(const SourceLinkAdapterIterator& other) const {
    return m_iterator - other.m_iterator;
  }

  BaseIterator m_iterator;
};

/// Delegate to unpack the surface associated with a source link
using SourceLinkSurfaceAccessor = Delegate<const Surface*(const SourceLink&)>;

}  // namespace Acts

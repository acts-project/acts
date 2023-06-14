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
#include "Acts/Utilities/TypeTraits.hpp"

#include <cassert>
#include <iostream>
#include <type_traits>
#include <utility>

namespace Acts {

namespace detail_sl {
template <typename T>
using geometry_id_t = decltype(std::declval<T>().geometryId());
}  // namespace detail_sl

class SourceLink final {
  using any_type = AnyBase<16>;

 public:
  /// Getter for the geometry identifier
  /// @return The GeometryIdentifier
  constexpr GeometryIdentifier geometryId() const { return m_geometryId; }

  SourceLink(const SourceLink& other) = default;
  SourceLink(SourceLink&& other) = default;
  SourceLink& operator=(const SourceLink& other) = default;
  SourceLink& operator=(SourceLink&& other) = default;

  /// Constructor from source link and explicit geometry id
  /// @tparam T The source link type
  /// @param id The geometry identifier
  /// @param upstream The upstream source link to store
  template <typename T, typename = std::enable_if_t<
                            !std::is_same_v<std::decay_t<T>, SourceLink>>>
  explicit SourceLink(GeometryIdentifier id, T&& upstream)
      : m_geometryId{id}, m_upstream{std::move(upstream)} {
    static_assert(!std::is_same_v<std::decay_t<T>, SourceLink>,
                  "Cannot wrap SourceLink in SourceLink");
  }

  /// Constructor from source link only, geometry identifier is determined
  /// automatically
  /// @tparam T The source link type
  /// @param upstream The upstream source link to store
  template <typename T, typename = std::enable_if_t<
                            Concepts::exists<detail_sl::geometry_id_t, T> &&
                            !std::is_same_v<std::decay_t<T>, SourceLink>>>
  explicit SourceLink(T&& upstream) : m_geometryId{upstream.geometryId()} {
    static_assert(
        std::is_same_v<detail_sl::geometry_id_t<T>, GeometryIdentifier>,
        "geometryId method does not return a geometry id type");

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
  GeometryIdentifier m_geometryId{};
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

}  // namespace Acts

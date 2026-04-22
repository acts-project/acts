// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/core/detail/data_context.hpp"
#include "detray/core/detail/tuple_container.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/utils/type_registry.hpp"
#include "detray/utils/type_traits.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <type_traits>

namespace detray {

/// @brief Wraps a vecmem enabled tuple and adds functionality to handle data
/// collections @tparam Ts.
///
/// @tparam An enum of type IDs that needs to match the value types of the Ts
/// pack.
/// @tparam context_t How to retrieve data according to e.g. conditions data
/// @tparam tuple_t The type of the underlying tuple container.
/// @tparam Ts the data collection types
template <typename ID = std::size_t, typename context_t = empty_context,
          template <typename...> class tuple_t = dtuple, typename... Ts>
class multi_store {
  /// Index range type that can be used with a typed index
  /// (i.e. leaves enough bits in the front unencoded for the type id)
  using index_range_t =
      detail::index_range<dindex, detail::sized_index_range,
                          std::uint_least32_t, 0x0fffff00, 0x000000ff>;

 public:
  using size_type = typename detail::first_t<Ts...>::size_type;
  using context_type = context_t;

  /// How to find and index a data collection in the store
  /// @{
  using ids = ID;
  using single_link = dtyped_index<ID, dindex>;
  using range_link = dtyped_index<ID, index_range_t>;
  /// @}

  /// Allow matching between IDs and collection value types
  /// @{
  using value_types = types::registry<ID, detail::get_value_t<Ts>...>;
  /// @}

  /// Underlying tuple container that can handle vecmem views
  using tuple_type = detail::tuple_container<tuple_t, Ts...>;
  /// Vecmem view types
  using view_type = typename tuple_type::view_type;
  using const_view_type = typename tuple_type::const_view_type;
  using buffer_type = typename tuple_type::buffer_type;

  /// Empty container
  constexpr multi_store() = default;
  /// Move constructor
  constexpr multi_store(multi_store &&) noexcept = default;

  // Delegate constructors to tuple container, which handles the memory

  /// Copy construct from element types
  constexpr explicit multi_store(const Ts &...args)
      : m_tuple_container(args...) {}

  /// Construct with a specific vecmem memory resource @param resource
  /// (host-side only)
  template <typename allocator_t = vecmem::memory_resource>
    requires(std::derived_from<allocator_t, std::pmr::memory_resource>)
  DETRAY_HOST explicit multi_store(allocator_t &resource)
      : m_tuple_container(resource) {}

  /// Copy Construct with a specific (vecmem) memory resource @param resource
  /// (host-side only)
  template <typename allocator_t = vecmem::memory_resource,
            typename T = tuple_t<Ts...>>
    requires(std::is_same_v<T, dtuple<Ts...>> &&
             std::derived_from<allocator_t, std::pmr::memory_resource>)
  DETRAY_HOST explicit multi_store(allocator_t &resource, const Ts &...args)
      : m_tuple_container(resource, args...) {}

  /// Construct from the container @param view . Mainly used device-side.
  template <concepts::device_view tuple_view_t>
  DETRAY_HOST_DEVICE explicit multi_store(tuple_view_t &view)
      : m_tuple_container(view) {}

  /// Move assignment operator
  multi_store &operator=(multi_store &&) noexcept = default;

  /// @returns a pointer to the underlying tuple container - const
  DETRAY_HOST_DEVICE
  constexpr auto data() const noexcept -> const tuple_type * {
    return &m_tuple_container;
  }

  /// @returns a pointer to the underlying tuple container - non-const
  DETRAY_HOST_DEVICE
  constexpr auto data() noexcept -> tuple_type * { return &m_tuple_container; }

  /// @returns the size of the underlying tuple
  DETRAY_HOST_DEVICE
  static constexpr auto n_collections() -> std::size_t { return sizeof...(Ts); }

  /// @returns the collections iterator at the start position.
  template <ID id = ID{0}>
  DETRAY_HOST_DEVICE constexpr auto begin(const context_type & /*ctx*/ = {}) {
    return detail::get<types::index_cast<value_types, id>>(m_tuple_container)
        .begin();
  }

  /// @returns the collections iterator sentinel.
  template <ID id = ID{0}>
  DETRAY_HOST_DEVICE constexpr auto end(const context_type & /*ctx*/ = {}) {
    return detail::get<types::index_cast<value_types, id>>(m_tuple_container)
        .end();
  }

  /// @returns a data collection by @tparam ID - const
  template <ID id>
  DETRAY_HOST_DEVICE constexpr decltype(auto) get(
      const context_type & /*ctx*/ = {}) const noexcept {
    return detail::get<types::index_cast<value_types, id>>(m_tuple_container);
  }

  /// @returns a data collection by @tparam ID - non-const
  template <ID id>
  DETRAY_HOST_DEVICE constexpr decltype(auto) get(
      const context_type & /*ctx*/ = {}) noexcept {
    return detail::get<types::index_cast<value_types, id>>(m_tuple_container);
  }

  /// @returns the size of a data collection by @tparam ID
  template <ID id>
  DETRAY_HOST_DEVICE constexpr auto size(
      const context_type & /*ctx*/ = {}) const noexcept -> dindex {
    return static_cast<dindex>(
        detail::get<types::index_cast<value_types, id>>(m_tuple_container)
            .size());
  }

  /// @returns the number of elements in all collections
  template <std::size_t current_idx = 0>
  DETRAY_HOST_DEVICE auto total_size(const context_type &ctx = {},
                                     dindex n = 0u) const noexcept -> dindex {
    n += size<types::id_cast<value_types, current_idx>>(ctx);

    if constexpr (current_idx < sizeof...(Ts) - 1) {
      return total_size<current_idx + 1>(ctx, n);
    }
    return n;
  }

  /// @returns true if the collection given by @tparam ID is empty
  template <ID id>
  DETRAY_HOST_DEVICE constexpr auto empty(
      const context_type & /*ctx*/ = {}) const noexcept -> bool {
    return detail::get<types::index_cast<value_types, id>>(m_tuple_container)
        .empty();
  }

  /// @returns true if every collection in container is empty.
  template <std::size_t current_idx = 0>
  DETRAY_HOST_DEVICE constexpr bool all_empty(const context_type &ctx = {},
                                              bool is_empty = true) const {
    is_empty &= empty<types::id_cast<value_types, current_idx>>(ctx);

    if constexpr (current_idx < sizeof...(Ts) - 1) {
      return all_empty<current_idx + 1>(ctx, is_empty);
    }
    return is_empty;
  }

  /// Removes and destructs all elements in a specific collection.
  template <ID id>
  DETRAY_HOST void clear(const context_type & /*ctx*/) {
    detail::get<types::index_cast<value_types, id>>(m_tuple_container).clear();
  }

  /// Removes and destructs all elements in the container.
  template <std::size_t current_idx = 0>
  DETRAY_HOST void clear_all(const context_type &ctx = {}) {
    clear<types::id_cast<value_types, current_idx>>(ctx);

    if constexpr (current_idx < sizeof...(Ts) - 1) {
      clear_all<current_idx + 1>(ctx);
    }
  }

  /// Reserve memory of size @param n for a collection given by @tparam id
  template <ID id>
  DETRAY_HOST void reserve(std::size_t n, const context_type & /*ctx*/) {
    detail::get<types::index_cast<value_types, id>>(m_tuple_container)
        .reserve(n);
  }

  /// Resize the underlying container to @param n for a collection given by
  /// @tparam id
  template <ID id>
  DETRAY_HOST void resize(std::size_t n, const context_type & /*ctx*/) {
    detail::get<types::index_cast<value_types, id>>(m_tuple_container)
        .resize(n);
  }

  /// Add a new element to a collection
  ///
  /// @tparam ID is the id of the collection
  /// @tparam T The element to be added to the collection
  ///
  /// @param args is the list of constructor arguments
  ///
  /// @note in general can throw an exception
  template <ID id, typename T>
  DETRAY_HOST constexpr auto push_back(
      const T &arg, const context_type & /*ctx*/ = {}) noexcept(false) -> void {
    auto &coll =
        detail::get<types::index_cast<value_types, id>>(m_tuple_container);
    return coll.push_back(arg);
  }

  /// Add a new element to a collection in place
  ///
  /// @tparam ID is the id of the collection
  /// @tparam Args are the types of the constructor arguments
  ///
  /// @param args is the list of constructor arguments
  ///
  /// @note in general can throw an exception
  template <ID id, typename... Args>
  DETRAY_HOST constexpr decltype(auto) emplace_back(
      const context_type & /*ctx*/ = {}, Args &&...args) noexcept(false) {
    auto &coll =
        detail::get<types::index_cast<value_types, id>>(m_tuple_container);
    return coll.emplace_back(std::forward<Args>(args)...);
  }

  /// Add a new collection - move/copy
  ///
  /// @tparam derived_collection_t is the type of the collection
  ///
  /// @param new_data is the new collection to be added
  ///
  /// @note in general can throw an exception
  template <typename derived_collection_t>
  DETRAY_HOST auto insert(derived_collection_t &&new_data,
                          const context_type & /*ctx*/ = {}) noexcept(false)
      -> void {
    using collection_t = std::remove_cvref_t<derived_collection_t>;

    static_assert((std::is_same_v<collection_t, Ts> || ...) == true,
                  "The type is not included in the parameter pack.");

    auto &coll = detail::get<collection_t>(m_tuple_container);

    coll.reserve(coll.size() +
                 std::forward<derived_collection_t>(new_data).size());
    coll.insert(coll.end(),
                std::forward<derived_collection_t>(new_data).begin(),
                std::forward<derived_collection_t>(new_data).end());
  }

  /// Append another store to the current one
  ///
  /// @tparam current_idx is the index to start unrolling
  ///
  /// @param other The other store
  ///
  /// @note in general can throw an exception
  template <std::size_t current_idx = 0>
  DETRAY_HOST void append(const multi_store &other,
                          const context_type &ctx = {}) noexcept(false) {
    auto &coll = other.template get<types::id_cast<value_types, current_idx>>();
    insert(coll, ctx);

    if constexpr (current_idx < sizeof...(Ts) - 1) {
      append<current_idx + 1>(other);
    }
  }

  /// Append another store to the current one - move
  ///
  /// @tparam current_idx is the index to start unrolling
  ///
  /// @param other The other store
  ///
  /// @note in general can throw an exception
  template <std::size_t current_idx = 0>
  DETRAY_HOST void append(multi_store &&other,
                          const context_type &ctx = {}) noexcept(false) {
    auto &coll = other.template get<types::id_cast<value_types, current_idx>>();
    insert(std::move(coll), ctx);

    if constexpr (current_idx < sizeof...(Ts) - 1) {
      append<current_idx + 1>(std::move(other));
    }
  }

  /// Calls a functor with a every data collection as parameter.
  ///
  /// @tparam functor_t functor that will be called on the store.
  /// @tparam Args argument types for the functor
  ///
  /// @param args additional functor arguments
  ///
  /// @return the functor output
  template <typename functor_t, typename... Args>
  DETRAY_HOST_DEVICE decltype(auto) apply(Args &&...args) const {
    return m_tuple_container.template apply<functor_t>(
        std::forward<Args>(args)...);
  }

  /// Calls a functor with a specific data collection (given by ID).
  ///
  /// @tparam functor_t functor that will be called on the group.
  /// @tparam Args argument types for the functor
  ///
  /// @param id the element id
  /// @param args additional functor arguments
  ///
  /// @return the functor output
  template <typename functor_t, typename... Args>
  DETRAY_HOST_DEVICE decltype(auto) visit(const ID id, Args &&...args) {
    return m_tuple_container.template visit<functor_t>(
        static_cast<std::size_t>(id), std::forward<Args>(args)...);
  }

  /// Calls a functor with a specific element of a data collection
  /// (given by a link).
  ///
  /// @tparam functor_t functor that will be called on the group.
  /// @tparam Args argument types for the functor
  ///
  /// @param id the element id
  /// @param args additional functor arguments
  ///
  /// @return the functor output
  template <typename functor_t, typename link_t, typename... Args>
    requires(std::same_as<link_t, single_link> ||
             std::same_as<link_t, range_link>)
  DETRAY_HOST_DEVICE decltype(auto) visit(const link_t link,
                                          Args &&...args) const {
    return m_tuple_container.template visit<functor_t>(
        static_cast<std::size_t>(detail::get<0>(link)), detail::get<1>(link),
        std::forward<Args>(args)...);
  }

  /// Print the types that are in the store
  DETRAY_HOST
  static constexpr void print() {
    types::print<types::list<detail::get_value_t<Ts>...>>();
  }

  /// @return the view on this tuple container - non-const
  DETRAY_HOST auto get_data() -> view_type {
    return m_tuple_container.get_data();
  }

  /// @return the view on this tuple container - const
  DETRAY_HOST auto get_data() const -> const_view_type {
    return m_tuple_container.get_data();
  }

 private:
  /// The underlying tuple container implementation
  tuple_type m_tuple_container;
};

/// Helper type for a data store that uses a single collection container
template <typename ID, typename context_t, template <typename...> class tuple_t,
          template <typename...> class container_t, typename... Ts>
using regular_multi_store =
    multi_store<ID, context_t, tuple_t, container_t<Ts>...>;

}  // namespace detray

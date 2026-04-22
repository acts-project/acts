// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/utils/concepts.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/type_list.hpp"

// System include(s)
#include <type_traits>
#include <utility>

namespace detray {

namespace types {

/// @brief match types with indices/identifiers and vice versa.
///
/// @tparam IDs enum that references the types (not used in base class)
/// @tparam registered_types the types that is mapped to the ID enum values
template <concepts::type_id ID, typename... registered_types>
class registry {
 public:
  /// Make the type IDs accessible
  using id = ID;
  /// Make the registered types accessible
  using type_list = detray::types::list<registered_types...>;

  /// Conventions for some basic info
  enum : std::size_t {
    n_types = sizeof...(registered_types),
    e_any = sizeof...(registered_types),
    e_unknown = sizeof...(registered_types) + 1,
  };

  /// @returns whether a given index/identifier can be mapped to a type.
  /// @{
  /// @param type_idx the type index to be checked
  DETRAY_HOST_DEVICE static constexpr bool is_valid(
      const std::size_t type_idx) {
    return type_idx < n_types;
  }

  /// @param type_id the type identifier to be checked
  template <typename id_t = id>
    requires(!std::convertible_to<id_t, std::size_t>)
  DETRAY_HOST_DEVICE static constexpr bool is_valid(const id_t type_id) {
    return static_cast<std::size_t>(type_id) < n_types;
  }
  /// @}
};

}  // namespace types

namespace concepts {

/// Type registry concepts
/// @{
template <typename R>
concept type_registry =
    type_list<typename R::type_list> && type_id<typename R::id>;
///@}

}  // namespace concepts

namespace types {

/// @brief Select a number of types from another registry and map their IDs to a
/// continuous range of new type indices.
///
/// @note Contrary to the original type registry, the type ID does not match the
/// position of the filtered type in the mapped registry anymore
///
/// @tparam registry_t the original type registry
/// @tparam type_selector_t how to select the new types from the original ones
template <concepts::type_registry registry_t, class type_selector_t>
class mapped_registry : public registry_t {
 public:
  /// Make the original type ids accessible
  using id = typename registry_t::id;
  /// Make the original types accessible
  using orig_types = typename registry_t::type_list;
  /// Make the filtered types accessible
  using type_list =
      decltype(detray::types::filtered_list<orig_types, type_selector_t>());

  /// Conventions for some basic info
  enum : std::size_t {
    n_types = detray::types::size<type_list>,
    e_any = detray::types::size<type_list>,
    e_unknown = detray::types::size<type_list> + 1u,
  };

  /// @returns whether a given index/identifier can be mapped to a filtered
  /// type.
  /// @{
  /// @param filtered_type_idx the filtered type index to be checked
  DETRAY_HOST_DEVICE static constexpr bool is_valid(
      const std::size_t filtered_type_idx) {
    return filtered_type_idx < n_types;
  }

  /// @param orig_id the original type identifier (ID) to be checked
  template <typename id_t = id>
    requires(!std::convertible_to<id_t, std::size_t>)
  DETRAY_HOST_DEVICE static constexpr bool is_valid(const id_t orig_id) {
    if (static_cast<std::size_t>(orig_id) >= index_map().size()) {
      return false;
    } else {
      return mapped_index(orig_id) < n_types;
    }
  }
  /// @}

  /// @returns the array that maps the original type positions to the new ones
  DETRAY_HOST_DEVICE
  static consteval auto index_map() {
    return detray::types::filtered_indices<orig_types, type_selector_t>(
        detray::types::list{});
  }

  /// @returns the filtered type index corresponding to the original type
  /// index
  DETRAY_HOST_DEVICE
  static constexpr std::size_t mapped_index(const std::size_t orig_idx) {
    return index_map()[orig_idx];
  }

  /// @returns the filtered type index corresponding to the original type ID
  DETRAY_HOST_DEVICE
  static constexpr std::size_t mapped_index(const id orig_id) {
    return index_map()[static_cast<std::size_t>(orig_id)];
  }

  /// @returns the original indexd given the new @param pos in the type list
  DETRAY_HOST_DEVICE
  static constexpr std::size_t original_index(const std::size_t pos) {
    std::size_t j{0u};
    for (std::size_t i : index_map()) {
      if (i == pos) {
        assert(registry_t::is_valid(j));
        return j;
      }
      j++;
    }
    return detray::detail::invalid_value<std::size_t>();
  }

  /// @returns the original id given the new @param pos in the type list
  DETRAY_HOST_DEVICE
  static constexpr id original_id(const std::size_t pos) {
    return static_cast<id>(original_index(pos));
  }
};

}  // namespace types

namespace concepts {

/// Mapped type registry concepts
/// @{
template <typename R>
concept mapped_type_registry =
    type_list<typename R::type_list> && type_list<typename R::orig_types> &&
    type_id<typename R::id>;
///@}

}  // namespace concepts

namespace types {

/// Specialization of the type list @c size traits
/// @see detray/utils/type_list.hpp
/// @{
template <concepts::type_id ID, typename... Ts>
struct get_size<registry<ID, Ts...>>
    : std::integral_constant<std::size_t, sizeof...(Ts)> {};

template <typename R, class S>
struct get_size<mapped_registry<R, S>>
    : std::integral_constant<std::size_t, mapped_registry<R, S>::n_types> {};
/// @}

/// Specialization of the type list @c contains trait
/// @see detray/utils/type_list.hpp
/// @{
template <typename T, std::size_t I, concepts::type_id ID, typename... Ts>
struct contains_impl<T, I, registry<ID, Ts...>>
    : public contains_impl<T, I, typename registry<ID, Ts...>::type_list> {};

template <typename T, std::size_t I, typename R, class S>
struct contains_impl<T, I, mapped_registry<R, S>> {
  // Figure out which type list to check
  // (use original list if 'T' is from there, otherwise try filtered list)
  using orig_list_t = typename mapped_registry<R, S>::orig_types;
  using filt_list_t = typename mapped_registry<R, S>::type_list;
  using type_list_t =
      std::conditional_t<contains<orig_list_t, T>, orig_list_t, filt_list_t>;

  using base_t = contains_impl<T, I, type_list_t>;

  static constexpr bool value = base_t::value;
  // The type is not contained in 'type_list'm the position will be invalid:
  // Do not pass that to the range-checked array of the filtered indices
  static constexpr std::size_t pos =
      (!base_t::value
           ? std::numeric_limits<std::size_t>::max()
           // Depending on which type list the type belongs to, map the index
           : (std::same_as<type_list_t, orig_list_t>
                  ? mapped_registry<R, S>::mapped_index(base_t::pos)
                  : base_t::pos));
};
/// @}

/// Specialization of the type list @c size trait
/// @see detray/utils/type_list.hpp
/// @{
template <std::size_t N, concepts::type_id ID, typename... Ts>
struct get_at<N, registry<ID, Ts...>> {
  static_assert(registry<ID, Ts...>::is_valid(N),
                "Index out of range in type registry 'at'");
  using type =
      typename get_at<N, typename registry<ID, Ts...>::type_list>::type;
};

template <std::size_t N, typename R, class S>
struct get_at<N, mapped_registry<R, S>> {
  static_assert(mapped_registry<R, S>::is_valid(N),
                "Index out of range in mapped type registry 'at'");
  using type =
      typename get_at<N, typename mapped_registry<R, S>::type_list>::type;
};
/// @}

/// Traits for the type registry
/// @{

/// Get the ID of a type
/// @{
template <typename = void, typename = void>
struct get_id {};

template <typename T, concepts::type_id ID, typename... Ts>
struct get_id<T, registry<ID, Ts...>> {
  static constexpr ID value =
      static_cast<ID>(position<typename registry<ID, Ts...>::type_list, T>);

  static_assert(contains<typename registry<ID, Ts...>::type_list, T>,
                "Unknown type in 'get id' for type registry");
  static_assert(registry<ID, Ts...>::is_valid(value),
                "Unmatched type ID in type registry 'id'");
};

template <typename T, typename R, class S>
struct get_id<T, mapped_registry<R, S>> {
  // Figure out which type list to check
  // (use original list if 'T' is from there, otherwise try filtered list)
  using mapped_reg_t = mapped_registry<R, S>;
  using orig_list_t = typename mapped_reg_t::orig_types;
  using filt_list_t = typename mapped_reg_t::type_list;
  using type_list_t =
      std::conditional_t<contains<orig_list_t, T>, orig_list_t, filt_list_t>;

  static_assert(std::same_as<type_list_t, orig_list_t> ||
                    (std::same_as<type_list_t, filt_list_t> &&
                     contains<filt_list_t, T>),
                "Unknown type in 'get id' for mapped type registry");

  using id_t = typename mapped_reg_t::id;
  static constexpr id_t value =
      contains<orig_list_t, T>
          ? static_cast<id_t>(position<orig_list_t, T>)
          : mapped_reg_t::original_id(position<filt_list_t, T>);
};

template <typename R, typename T>
inline constexpr R::id id{get_id<T, R>::value};
/// @}

/// Get the type corresponding to an ID
/// @{
template <std::size_t I, typename = void>
struct get_type {};

template <std::size_t I, typename ID, typename... Ts>
struct get_type<I, registry<ID, Ts...>> {
  static_assert(registry<ID, Ts...>::is_valid(I),
                "Unmatched type ID in type registry 'get'");

  using type = at<registry<ID, Ts...>, I>;
};

template <std::size_t I, typename R, class S>
struct get_type<I, mapped_registry<R, S>> {
  static_assert(
      mapped_registry<R, S>::is_valid(mapped_registry<R, S>::mapped_index(I)),
      "Unmatched type ID in mapped type registry 'get'");

  using type =
      at<mapped_registry<R, S>, mapped_registry<R, S>::mapped_index(I)>;
};

template <typename R, auto I>
using get = typename get_type<static_cast<std::size_t>(I), R>::type;
/// @}

/// Get the type corresponding to an ID
/// @{
template <std::size_t I, typename = void>
struct cast_impl_id {};

template <std::size_t I, typename ID, typename... Ts>
struct cast_impl_id<I, registry<ID, Ts...>> {
  static_assert(registry<ID, Ts...>::is_valid(I),
                "Cannot cast to ID: Invalid index for type registry");

  static constexpr ID value = static_cast<ID>(I);
};

template <std::size_t I, typename R, class S>
struct cast_impl_id<I, mapped_registry<R, S>> {
  static_assert(
      mapped_registry<R, S>::is_valid(mapped_registry<R, S>::mapped_index(I)),
      "Cannot cast to ID: Invalid index for ,apped type registry");

  using ID = typename mapped_registry<R, S>::id;
  static constexpr ID value = static_cast<ID>(I);
};

template <std::size_t I, typename = void>
struct cast_impl_idx {};

template <std::size_t I, typename ID, typename... Ts>
struct cast_impl_idx<I, registry<ID, Ts...>> {
  static_assert(registry<ID, Ts...>::is_valid(I),
                "Cannot cast to index: Invalid index for type registry");

  static constexpr std::size_t value = I;
};

template <std::size_t I, typename R, class S>
struct cast_impl_idx<I, mapped_registry<R, S>> {
  static_assert(
      mapped_registry<R, S>::is_valid(mapped_registry<R, S>::mapped_index(I)),
      "Cannot cast to index: Invalid index for ,apped type registry");

  static constexpr std::size_t value = mapped_registry<R, S>::mapped_index(I);
};

template <typename R, std::size_t I>
inline constexpr R::id id_cast = cast_impl_id<I, R>::value;

template <typename R, auto I>
inline constexpr std::size_t index_cast =
    cast_impl_idx<static_cast<std::size_t>(I), R>::value;
/// @}

/// @}

/// Variadic unrolling of the tuple that calls a functor on the element that
/// corresponds to @param idx.
///
/// @tparam functor_t functor that will be called on the element.
/// @tparam Args argument types for the functor
/// @tparam current_idx Current index into the container tuple. Is converted
///         to an id_t and tested against the given id.
/// @tparam remaining_idcs te remaining tuple indices to be tested.
///
/// @see https://godbolt.org/z/qd6xns7KG
template <concepts::type_registry registry_t, typename functor_t,
          typename... Args, std::size_t current_idx = 0u,
          std::size_t... remaining_idcs>
  requires(std::invocable<functor_t,
                          const types::front<typename registry_t::type_list>&,
                          Args...>)
DETRAY_HOST_DEVICE constexpr decltype(auto) visit_helper(
    const std::size_t idx,
    std::index_sequence<current_idx, remaining_idcs...> /*seq*/,
    Args&&... args) {
  // Check if the current tuple index is matched to the target index
  if constexpr (concepts::mapped_type_registry<registry_t>) {
    if (registry_t::mapped_index(idx) == current_idx) {
      return functor_t{}(types::at<registry_t, current_idx>{},
                         std::forward<Args>(args)...);
    }
  } else {
    if (idx == current_idx) {
      return functor_t{}(types::at<registry_t, current_idx>{},
                         std::forward<Args>(args)...);
    }
  }

  // Check the next index
  using return_t = std::invoke_result_t<
      functor_t, const types::front<typename registry_t::type_list>&, Args...>;

  if constexpr (sizeof...(remaining_idcs) >= 1u) {
    using next_elem_t = types::at<registry_t, current_idx + 1u>;
    using next_return_t =
        std::invoke_result_t<functor_t, const next_elem_t&, Args...>;
    static_assert(
        std::same_as<return_t, next_return_t>,
        "Functor return type must be the same for all elements of the "
        "tuple.");

    return visit_helper<registry_t, functor_t>(
        idx, std::index_sequence<remaining_idcs...>{},
        std::forward<Args>(args)...);
  } else if constexpr (std::is_same_v<return_t, void>) {
    return;
  } else {
    return return_t{};
  }
}

/// Visits a tuple element according to its @param idx and calls
/// @tparam functor_t with the arguments @param args on it.
///
/// @returns the functor result (this is necessarily always of the same
/// type, regardless the input tuple element type).
template <concepts::type_registry registry_t, typename functor_t,
          typename... Args>
DETRAY_HOST_DEVICE constexpr decltype(auto) visit(const std::size_t idx,
                                                  Args&&... args) {
  return visit_helper<registry_t, functor_t>(
      idx, std::make_index_sequence<registry_t::n_types>{},
      std::forward<Args>(args)...);
}

/// Visits a tuple element according to its type_id @param id and calls
/// @tparam functor_t with the arguments @param args on it.
///
/// @returns the functor result (this is necessarily always of the same
/// type, regardless the input tuple element type).
template <concepts::type_registry registry_t, typename functor_t,
          typename... Args>
DETRAY_HOST_DEVICE constexpr decltype(auto) visit(
    const typename registry_t::id type_id, Args&&... args) {
  assert(registry_t::is_valid(type_id));

  return visit<registry_t, functor_t>(static_cast<std::size_t>(type_id),
                                      std::forward<Args>(args)...);
}

}  // namespace types

}  // namespace detray

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/detail/compact_device_vector.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/tuple.hpp"

// Vecmem include(s)
#include <vecmem/containers/jagged_vector.hpp>
#include <vecmem/containers/vector.hpp>

// System include(s)
#include <array>
#include <map>
#include <type_traits>
#include <vector>

namespace detray {

template <typename value_t, std::size_t kDIM>
/// Fixed-size array container.
using darray = std::array<value_t, kDIM>;

template <typename value_t>
/// Vector container backed by vecmem.
using dvector = vecmem::vector<value_t>;

template <typename value_t>
/// Jagged vector container backed by vecmem.
using djagged_vector = vecmem::jagged_vector<value_t>;

template <typename key_t, typename value_t>
/// Ordered key-value map container.
using dmap = std::map<key_t, value_t>;

template <class... types>
/// Tuple type used by detray storage and dispatch.
using dtuple = detray::tuple<types...>;

/// @brief Bundle container type definitions
template <template <typename...> class vector_t = dvector>
struct container_types {
  template <typename T>
  /// Vector template selected for this container policy.
  using vector_type = vector_t<T>;

  template <class... T>
  /// Tuple template used to combine component types.
  using tuple_type = dtuple<T...>;

  template <typename T, std::size_t kDIM>
  /// Fixed-size array template used for inline collections.
  using array_type = darray<T, kDIM>;
};

/// @brief Bundle container type definitions
template <template <typename...> class vector_t = dvector>
struct const_container_types {
  template <typename T>
  /// Vector template selected for this container policy.
  using vector_type = vector_t<std::add_const_t<T>>;

  template <class... T>
  /// Tuple template used to combine component types.
  using tuple_type = dtuple<std::add_const_t<T>...>;

  template <typename T, std::size_t kDIM>
  /// Fixed-size array template used for inline collections.
  using array_type = darray<std::add_const_t<T>, kDIM>;
};

/// Defining some common types
/// Default host container policy based on vecmem vectors.
using host_container_types = container_types<>;
/// Const-qualified host container policy based on vecmem vectors.
using const_host_container_types = const_container_types<>;

namespace detail {

// make std::get available in detray detail namespace, where also the thrust and
// index specific overloads live.
using std::get;

/// Trait class to figure out if a given type has a @c reserve(...) function
template <typename T>
struct has_reserve {
 private:
  /// Function returning @c std::true_type for types that do have a @c
  /// reserve(...) function
  template <typename C>
  static constexpr auto check(C*) ->
      typename std::is_void<decltype(std::declval<C>().reserve(
          std::declval<typename C::size_type>()))>::type;

  /// Function returning @c std::false_type for types that fair the previous
  /// function
  template <typename>
  static constexpr std::false_type check(...);

  /// Declare the value type of this trait class
  using type = decltype(check<T>(nullptr));

 public:
  /// Value of the check
  static constexpr bool value = type::value;
};

/// @name Functions calling or not calling reserve(...) based on whether it's
/// available
/// @{
template <typename T>
  requires has_reserve<T>::value
DETRAY_HOST_DEVICE void call_reserve(T& obj, std::size_t newsize) {
  obj.reserve(newsize);
}

template <typename T>
  requires(!has_reserve<T>::value)
DETRAY_HOST_DEVICE void call_reserve(T& /*obj*/, std::size_t /*newsize*/) {
  /*Not defined*/
}
/// @}

}  // namespace detail

}  // namespace detray

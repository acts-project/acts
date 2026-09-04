// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/tuple.hpp"
#include "detray/utils/tuple_helpers.hpp"

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>

// System include(s)
#include <concepts>

namespace detray {

/// Container types used in device code
using device_container_types =
    container_types<compact_device_vector, vecmem::jagged_device_vector>;

/// Specialized view for @c vecmem::vector containers
template <typename T>
using dvector_view = vecmem::data::vector_view<T>;

/// Specialized view for @c vecmem::jagged_vector containers
template <typename T>
using djagged_vector_view = vecmem::data::jagged_vector_view<T>;

namespace detail {
// Views for types that aggregate containers/other viewable types

/// Empty view type for inheritance template resolution
struct dbase_view {};

/// Helper trait to determine if a type can be interpreted as a (composite)
/// vecemem view
/// @{
template <typename T>
struct is_device_view : public std::false_type {};

/// Specialization of 'is view' for @c vecmem::data::vector_view containers
template <typename T>
struct is_device_view<vecmem::data::vector_view<T>> : public std::true_type {};

/// Specialization of 'is view' for constant @c vecmem::data::vector_view
/// containers
template <typename T>
struct is_device_view<const vecmem::data::vector_view<T>>
    : public std::true_type {};

/// Specialization of 'is view' for @c vecmem::data::jagged_vector_view
/// containers
template <typename T>
struct is_device_view<vecmem::data::jagged_vector_view<T>>
    : public std::true_type {};

/// Specialization of 'is view' for constant @c vecmem::data::jagged_vector_view
/// containers
template <typename T>
struct is_device_view<const vecmem::data::jagged_vector_view<T>>
    : public std::true_type {};

template <typename T>
inline constexpr bool is_device_view_v = is_device_view<T>::value;
/// @}
}  // namespace detail

namespace concepts {

template <typename T>
concept device_view =
    std::derived_from<std::remove_cvref_t<T>, detray::detail::dbase_view> ||
    detail::is_device_view_v<T>;

}  // namespace concepts

namespace detail {
/// Helper trait to check whether a type has a vecmem view defined
/// @{
template <class T>
struct get_view : public std::false_type {
  using type = void;
};

/// Get view of detray type
template <class T>
  requires concepts::device_view<typename T::view_type>
struct get_view<T> : public std::true_type {
  using type =
      std::conditional_t<std::is_const_v<T>, typename T::const_view_type,
                         typename T::view_type>;
};

/// Specialization of the view getter for @c vecmem::vector
template <typename T>
struct get_view<vecmem::vector<T>> : public std::true_type {
  using type = dvector_view<T>;
};

/// Specialization of the view getter for @c vecmem::vector - const
template <typename T>
struct get_view<const vecmem::vector<T>> : public std::true_type {
  using type = dvector_view<const T>;
};

/// Specialization of the view getter for @c vecmem::device_vector
template <typename T>
struct get_view<vecmem::device_vector<T>> : public std::true_type {
  using type = dvector_view<T>;
};

/// Specialization of the view getter for @c vecmem::device_vector - const
template <typename T>
struct get_view<const vecmem::device_vector<T>> : public std::true_type {
  using type = dvector_view<const T>;
};

/// Specialization of the view getter for @c vecmem::device_vector
template <typename T>
struct get_view<compact_device_vector<T>> : public std::true_type {
  using type = dvector_view<T>;
};

/// Specialization of the view getter for @c vecmem::device_vector - const
template <typename T>
struct get_view<const compact_device_vector<T>> : public std::true_type {
  using type = dvector_view<const T>;
};

/// Specialization of the view getter for @c vecmem::jagged_vector
template <typename T>
struct get_view<vecmem::jagged_vector<T>> : public std::true_type {
  using type = djagged_vector_view<T>;
};

/// Specialization of the view getter for @c vecmem::jagged_vector - const
template <typename T>
struct get_view<const vecmem::jagged_vector<T>> : public std::true_type {
  using type = djagged_vector_view<const T>;
};

/// Specialization of the view getter for @c vecmem::jagged_device_vector
template <typename T>
struct get_view<vecmem::jagged_device_vector<T>> : public std::true_type {
  using type = djagged_vector_view<T>;
};

/// Specialization of the view getter for @c vecmem::jagged_device_vector
/// - const
template <typename T>
struct get_view<const vecmem::jagged_device_vector<T>> : public std::true_type {
  using type = djagged_vector_view<const T>;
};

template <class T>
using get_view_t = typename get_view<T>::type;
/// @}
}  // namespace detail

namespace concepts {

template <class T>
concept viewable = concepts::device_view<detail::get_view_t<T>>;

}  // namespace concepts

/// @brief General view type that aggregates vecmem based view implementations.
///
/// This is for detray types that hold multiple members that all define custom
/// view types of their own. The 'sub'-views are being aggregated in this helper
/// and are extracted in the types constructor and then handed down to the
/// member constructors.
template <concepts::device_view... view_ts>
struct dmulti_view : public detail::dbase_view {
  dtuple<view_ts...> m_view;

  dmulti_view() = default;

  /// Tie multiple views together
  DETRAY_HOST_DEVICE
  explicit dmulti_view(view_ts&&... views) : m_view(std::move(views)...) {}

  /// Tie multiple views together
  DETRAY_HOST_DEVICE
  explicit dmulti_view(view_ts&... views) : m_view(views...) {}
};

/// @brief Detray version of 'get_data' - non-const
///
/// It is available to the generic containers before the value types are known,
/// thus enabling 'duck typing'.
///
/// @note This does not pick up the vecmem types.
template <concepts::viewable T>
typename T::view_type get_data(T& viewable) {
  return viewable.get_data();
}

/// @brief Detray version of 'get_data' - const
///
/// It is available to the generic containers before the value types are known,
/// thus enabling 'duck typing'.
///
/// @note This does not pick up the vecmem types.
template <concepts::viewable T>
typename T::const_view_type get_data(const T& viewable) {
  return viewable.get_data();
}

/// Get the vecmem view type of a vector
template <typename T, typename A>
dvector_view<T> get_data(std::vector<T, A>& vec) {
  return vecmem::get_data(vec);
}

/// Get the vecmem view type of a vector - const
template <typename T, typename A>
dvector_view<const T> get_data(const std::vector<T, A>& vec) {
  return vecmem::get_data(vec);
}

}  // namespace detray

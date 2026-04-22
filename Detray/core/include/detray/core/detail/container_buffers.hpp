// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/tuple.hpp"
#include "detray/utils/tuple_helpers.hpp"

// Vecmem include(s)
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

// System include(s)
#include <concepts>
#include <type_traits>

namespace detray {

/// Flag how a copy should be done, only needed until copy util is ready
enum class copy {
  sync = 0,
  async = 1,
};

namespace detail {
// Buffer wrappers for types that aggregate containers/other bufferable types

/// Empty buffer type for inheritance template resolution
struct dbase_buffer {};

/// Helper trait to determine if a type can be interpreted as a (composite)
/// vecemem buffer. This is the case if it inherits from @c dbase_buffer or if
/// it matches one of the vecmem specializations
/// @{
template <typename T>
struct is_buffer : public std::false_type {};

/// Specialization of @c is_buffer for types that inherit from @c dbase_buffer
template <typename T>
  requires std::is_base_of_v<detray::detail::dbase_buffer,
                             std::remove_cvref_t<T>>
struct is_buffer<T> : public std::true_type {};

/// Specialization of @c is_buffer for @c vecmem::data::vector_buffer containers
template <typename T>
struct is_buffer<vecmem::data::vector_buffer<T>> : public std::true_type {};

/// Specialization of @c is_buffer for constant @c vecmem::data::vector_buffer
/// containers
template <typename T>
struct is_buffer<const vecmem::data::vector_buffer<T>> : public std::true_type {
};

template <typename T>
inline constexpr bool is_buffer_v = is_buffer<T>::value;
/// @}

}  // namespace detail

namespace concepts {

template <typename T>
concept device_buffer =
    std::derived_from<std::remove_cvref_t<T>, detray::detail::dbase_buffer> ||
    detray::detail::is_buffer_v<T>;

}  // namespace concepts

namespace detail {
/// Helper trait to check whether a type has a [vecmem] buffer type defined
/// @{
template <class T>
struct get_buffer {};

template <class T>
  requires detray::detail::is_buffer_v<
      typename std::remove_cvref_t<T>::buffer_type>
struct get_buffer<T> {
  using type = typename T::buffer_type;
};

/// Specialization of the buffer getter for @c vecmem::vector
template <typename T>
struct get_buffer<vecmem::vector<T>> {
  using type = vecmem::data::vector_buffer<T>;
};

/// Specialization of the buffer getter for @c vecmem::vector - const
template <typename T>
struct get_buffer<const vecmem::vector<T>> {
  using type = vecmem::data::vector_buffer<const T>;
};

/// Specialization of the buffer getter for @c vecmem::device_vector
template <typename T>
struct get_buffer<vecmem::device_vector<T>> {
  using type = vecmem::data::vector_buffer<T>;
};

/// Specialization of the buffer getter for @c vecmem::device_vector - const
template <typename T>
struct get_buffer<const vecmem::device_vector<T>> {
  using type = vecmem::data::vector_buffer<const T>;
};

/// Specialization of the buffer getter for @c vecmem::device_vector
template <typename T>
struct get_buffer<compact_device_vector<T>> {
  using type = vecmem::data::vector_buffer<T>;
};

/// Specialization of the buffer getter for @c vecmem::device_vector - const
template <typename T>
struct get_buffer<const compact_device_vector<T>> {
  using type = vecmem::data::vector_buffer<const T>;
};

template <class T>
using get_buffer_t = typename get_buffer<T>::type;
/// @}

}  // namespace detail

namespace concepts {

template <class T>
concept bufferable = concepts::device_buffer<detray::detail::get_buffer_t<T>>;

}  // namespace concepts
/// Specialized buffer for @c vecmem::vector containers
template <typename T>
using dvector_buffer = vecmem::data::vector_buffer<T>;

/// @brief General buffer type that hierarchically aggregates vecmem based
/// buffer implementations.
///
/// This is for detray classes that hold multiple members that all define custom
/// buffer types of their own. The 'sub'-buffers are being aggregated in this
/// helper and are extracted again in the @c get_data call when building the
/// buffer views. Since they are packaged into the @c dmulti_buffer_helper
/// according to the class member hierarchy, transcribing them into a view type
/// results back in a view type that is compatible with what needs to be passed
/// to the constructors of the original class.
///
/// class -> view -> buffer -> view (passed to kernel) -> device-side class
template <concepts::device_buffer... buffer_ts>
struct dmulti_buffer : public detail::dbase_buffer {
  dtuple<buffer_ts...> m_buffer;

  dmulti_buffer() = default;

  /// Tie multiple buffers together
  DETRAY_HOST
  explicit dmulti_buffer(buffer_ts&&... buffers)
      : m_buffer(std::move(buffers)...) {}
};

/// @brief Get the buffer representation of a vecmem vector - non-const
template <class T>
auto get_buffer(const dvector_view<T>& vec_view, vecmem::memory_resource& mr,
                vecmem::copy& cpy, detray::copy cpy_type = detray::copy::sync,
                vecmem::data::buffer_type buff_type =
                    vecmem::data::buffer_type::fixed_size) {
  // In case the view references a const object, return a non-const buffer
  using ret_buffer_t = dvector_buffer<std::remove_cv_t<T>>;

  ret_buffer_t buff{vec_view.size(), mr, buff_type};

  // TODO: Move this to detray copy util, which bundles vecmem copy object and
  // stream handle and gets this switch case right automatically
  if (cpy_type == detray::copy::async) {
    cpy(vec_view, buff)->ignore();
  } else {
    cpy(vec_view, buff)->wait();
  }

  return buff;
}

/// @brief Recursively get the buffer representation of a composite view
///
/// @note This does not pick up the vecmem types.
template <typename... Ts>
auto get_buffer(
    const dmulti_view<Ts...>&, vecmem::memory_resource&, vecmem::copy&,
    detray::copy = detray::copy::sync,
    vecmem::data::buffer_type =
        vecmem::data::buffer_type::fixed_size);  // Forward declaration

/// @brief Recursively get the buffer representation of a composite view
///
/// Unwraps the view type at compile time and calls @c get_buffer on every view.
/// Then returns the resulting buffer objects and packages them recursively
/// into a @c multi_buffer .
///
/// @note This does not pick up the vecmem types.
template <concepts::device_view... Ts, std::size_t... I>
auto get_buffer(const dmulti_view<Ts...>& data_view,
                vecmem::memory_resource& mr, vecmem::copy& cpy,
                std::index_sequence<I...> /*seq*/,
                detray::copy cpy_type = detray::copy::sync,
                vecmem::data::buffer_type buff_type =
                    vecmem::data::buffer_type::fixed_size) {
  // Evaluate recursive buffer type
  // (e.g. dmulti_view<..., dmulti_view<dvector_view<T>, ...>, ...>
  //       => dmulti_buffer<..., dmulti_buffer<dvector_buffer<T>, ...>, ...>)
  using result_buffer_t = dmulti_buffer<decltype(detray::get_buffer(
      detail::get<I>(data_view.m_view), mr, cpy, cpy_type, buff_type))...>;

  return result_buffer_t{std::move(detray::get_buffer(
      detail::get<I>(data_view.m_view), mr, cpy, cpy_type, buff_type))...};
}

template <typename... Ts>
auto get_buffer(const dmulti_view<Ts...>& data_view,
                vecmem::memory_resource& mr, vecmem::copy& cpy,
                detray::copy cpy_type, vecmem::data::buffer_type buff_type) {
  return detray::get_buffer(
      data_view, mr, cpy,
      std::make_index_sequence<
          detail::tuple_size_v<decltype(data_view.m_view)>>{},
      cpy_type, buff_type);
}

/// @brief Get the buffer representation of a composite object - non-const
template <concepts::bufferable T>
typename T::buffer_type get_buffer(T& bufferable, vecmem::memory_resource& mr,
                                   vecmem::copy& cpy,
                                   detray::copy cpy_type = detray::copy::sync,
                                   vecmem::data::buffer_type buff_type =
                                       vecmem::data::buffer_type::fixed_size) {
  return detray::get_buffer(bufferable.get_data(), mr, cpy, cpy_type,
                            buff_type);
}

/// @brief Get the buffer representation of a vecmem vector - non-const
template <class T>
dvector_buffer<std::remove_cv_t<T>> get_buffer(
    dvector<T>& vec, vecmem::memory_resource& mr, vecmem::copy& cpy,
    detray::copy cpy_type = detray::copy::sync,
    vecmem::data::buffer_type buff_type = vecmem::data::buffer_type::fixed_size
    /*, stream*/) {
  return detray::get_buffer(detray::get_data(vec), mr, cpy, cpy_type,
                            buff_type);
}

/// @brief Get the buffer representation of a vecmem vector - non-const
template <class T>
dvector_buffer<std::remove_cv_t<T>> get_buffer(
    const dvector<T>& vec, vecmem::memory_resource& mr, vecmem::copy& cpy,
    detray::copy cpy_type = detray::copy::sync,
    vecmem::data::buffer_type buff_type = vecmem::data::buffer_type::fixed_size
    /*, stream*/) {
  return detray::get_buffer(detray::get_data(vec), mr, cpy, cpy_type,
                            buff_type);
}

/// Get the vecmem view type of a vector buffer
///
/// @note this is needed, because the corresponding vecmem::get_data() function
///       would return a reference to a view.
template <class T>
dvector_view<T> get_data(dvector_buffer<T>& buff) {
  return vecmem::get_data(buff);
}

/// Get the vecmem view type of a vector buffer - const
template <class T>
dvector_view<const T> get_data(const dvector_buffer<T>& buff) {
  return vecmem::get_data(buff);
}

/// @brief Get the view ( @c dmulti_view ) of a @c dmulti_buffer
template <concepts::device_buffer... Ts>
auto get_data(dmulti_buffer<Ts...>& multi_buff);  // Forward declaration

/// @brief Get the view ( @c dmulti_view ) of a @c dmulti_buffer - const
template <concepts::device_buffer... Ts>
auto get_data(const dmulti_buffer<Ts...>& multi_buff);  // Forward declaration

/// @brief Unroll the composite buffer type
///
/// Unwraps the buffer type at compile time and calls @c get_data on every
/// buffer.
/// Then returns the resulting view objects and packages them into a
/// @c mutli_view to be ultimately passed on to the class constructor.
///
/// @note This does not pick up the vecmem types.
template <concepts::device_buffer... Ts, std::size_t... I>
auto get_data(dmulti_buffer<Ts...>& multi_buff, std::index_sequence<I...>) {
  // Evaluate recursive view type (reverse of 'get_buffer(dmulti_view)')
  using result_view_t =
      dmulti_view<decltype(detray::get_data(std::declval<Ts&>()))...>;

  return result_view_t{
      detray::get_data(detail::get<I>(multi_buff.m_buffer))...};
}

/// @brief Unroll the composite buffer type - const
template <concepts::device_buffer... Ts, std::size_t... I>
auto get_data(const dmulti_buffer<Ts...>& multi_buff,
              std::index_sequence<I...>) {
  // Evaluate recursive view type (reverse of 'get_buffer(dmulti_view)')
  using result_view_t =
      dmulti_view<decltype(detray::get_data(std::declval<const Ts&>()))...>;

  return result_view_t{
      detray::get_data(detail::get<I>(multi_buff.m_buffer))...};
}

template <concepts::device_buffer... Ts>
auto get_data(dmulti_buffer<Ts...>& multi_buff) {
  return detray::get_data(multi_buff,
                          std::make_index_sequence<sizeof...(Ts)>{});
}

template <concepts::device_buffer... Ts>
auto get_data(const dmulti_buffer<Ts...>& multi_buff) {
  return detray::get_data(multi_buff,
                          std::make_index_sequence<sizeof...(Ts)>{});
}

}  // namespace detray

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_type_list.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/geometry/move_only_any.hpp"

// Detray include(s).
#include <any>

namespace traccc {

class detector_buffer {
 public:
  /// @tparam T either a detray metadata or detector type (checked by the trait)
  template <typename T>
  void set(typename detray::detector_traits<T>::buffer&& obj) {
    m_obj.set<typename detray::detector_traits<T>::buffer>(std::move(obj));
  }

  /// @tparam T either a detray metadata or detector type (checked by the trait)
  template <typename T>
  bool is() const {
    return (type() == typeid(typename detray::detector_traits<T>::buffer));
  }

  const std::type_info& type() const { return m_obj.type(); }

  /// @tparam T either a detray metadata or detector type (checked by the trait)
  template <typename T>
  const typename detray::detector_traits<T>::buffer& as() const {
    return m_obj.as<typename detray::detector_traits<T>::buffer>();
  }

  /// @tparam T either a detray metadata or detector type (checked by the trait)
  template <typename T>
  typename detray::detector_traits<T>::view as_view() const {
    return detray::get_data(as<typename detray::detector_traits<T>::view>());
  }

 private:
  move_only_any m_obj;
};  // class detector_buffer

/// @brief Helper function for `detector_buffer_visitor`
template <typename callable_t,
          detray::concepts::detector_traits detector_traits_t,
          detray::concepts::detector_traits... detector_traits_ts>
auto detector_buffer_visitor_helper(
    const detector_buffer& detector_buffer, callable_t&& callable,
    std::tuple<detector_traits_t, detector_traits_ts...>*) {
  using detector_t = typename detector_traits_t::host;

  if (detector_buffer.is<detector_t>()) {
    return callable.template operator()<detector_t>(
        detector_buffer.as_view<detector_t>());
  } else {
    if constexpr (sizeof...(detector_traits_ts) > 0) {
      return detector_buffer_visitor_helper(
          detector_buffer, std::forward<callable_t>(callable),
          static_cast<std::tuple<detector_traits_ts...>*>(nullptr));
    } else {
      std::stringstream exception_message;

      exception_message << "Invalid detector type ("
                        << detector_buffer.type().name()
                        << ") received, but this type is not supported"
                        << std::endl;

      throw std::invalid_argument(exception_message.str());
    }
  }
}

/// @brief Visitor for polymorphic detector buffer types
///
/// This function takes a list of supported detector trait types and checks
/// if the provided field is one of them. If it is, it will call the provided
/// callable on a view of it and otherwise it will throw an exception.
template <typename detector_buffer_list_t, typename callable_t>
auto detector_buffer_visitor(const detector_buffer& detector_buffer,
                             callable_t&& callable) {
  return detector_buffer_visitor_helper(
      detector_buffer, std::forward<callable_t>(callable),
      static_cast<detector_buffer_list_t*>(nullptr));
}

// TODO: Docs
inline detector_buffer buffer_from_host_detector(const host_detector& det,
                                                 vecmem::memory_resource& mr,
                                                 vecmem::copy& copy) {
  return host_detector_visitor<traccc::detector_type_list>(
      det, [&mr, &copy]<detray::concepts::detector_traits detector_traits_t>(
               const typename detector_traits_t::host& detector) {
        traccc::detector_buffer rv;
        rv.set<typename detector_traits_t::host>(
            detray::get_buffer(detector, mr, copy));
        return rv;
      });
}

}  // namespace traccc

/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_type_list.hpp"
#include "traccc/geometry/move_only_any.hpp"

// Detray include(s).
#include <any>

namespace traccc {

/// Typeless, owning, host detector object
class host_detector {
 public:
  /// @tparam T either a detray metadata or detector type (checked by the trait)
  /// @{
  template <typename T>
  void set(detray::detector_host_t<T>&& obj) {
    m_obj.set<detray::detector_host_t<T>>(std::move(obj));
  }

  template <typename T>
  bool is() const {
    return (type() == typeid(detray::detector_host_t<T>));
  }

  const std::type_info& type() const { return m_obj.type(); }

  template <typename T>
  const detray::detector_host_t<T>& as() const {
    return m_obj.as<detray::detector_host_t<T>>();
  }
  /// @}

 private:
  move_only_any m_obj;
};

/// @brief Helper function for `host_detector_visitor`
template <typename callable_t, detray::concepts::detector detector_t,
          detray::concepts::detector... detector_ts>
auto host_detector_visitor_helper(const host_detector& host_detector,
                                  callable_t&& callable,
                                  std::tuple<detector_t, detector_ts...>*) {
  if (host_detector.is<detector_t>()) {
    return callable.template operator()<detector_t>(
        host_detector.as<detector_t>());
  } else {
    if constexpr (sizeof...(detector_ts) > 0) {
      return host_detector_visitor_helper(
          host_detector, std::forward<callable_t>(callable),
          static_cast<std::tuple<detector_ts...>*>(nullptr));
    } else {
      std::stringstream exception_message;

      exception_message << "Invalid detector type ("
                        << host_detector.type().name()
                        << ") received, but this type is not supported"
                        << std::endl;

      throw std::invalid_argument(exception_message.str());
    }
  }
}

/// @brief Visitor for polymorphic host detector types
///
/// This function takes a list of supported detector trait types and checks
/// if the provided field is one of them. If it is, it will call the provided
/// callable on a view of it and otherwise it will throw an exception.
template <typename detector_list_t, typename callable_t>
auto host_detector_visitor(const host_detector& host_detector,
                           callable_t&& callable) {
  return host_detector_visitor_helper(host_detector,
                                      std::forward<callable_t>(callable),
                                      static_cast<detector_list_t*>(nullptr));
}
}  // namespace traccc

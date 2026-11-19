// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "detray/core/host_detector.hpp"

namespace detray {

/// Owning, move-only detector buffer. Its memory resource must outlive it.
/// Views borrow the buffer; asynchronous work must finish before destruction.
class detector_buffer {
 public:
  template <typename T>
  void set(detector_buffer_t<T>&& object) {
    m_obj.set<detector_buffer_t<T>>(std::move(object));
  }
  template <typename T>
  bool is() const noexcept {
    return m_obj.is<detector_buffer_t<T>>();
  }
  template <typename T>
  const detector_buffer_t<T>& as() const {
    return m_obj.as<detector_buffer_t<T>>();
  }
  template <typename T>
  detector_view_t<T> as_view() const {
    return get_data(as<T>());
  }
  const std::type_info& type() const { return m_obj.type(); }
  bool has_value() const noexcept { return m_obj.has_value(); }
  void reset() noexcept { m_obj.reset(); }

 private:
  move_only_any m_obj;
};

/// Host dispatch to a concrete buffer view. The callable must not retain a view
/// beyond the lifetime of the buffer or its allocator.
template <typename detector_list_t, typename callable_t>
decltype(auto) detector_buffer_visitor(const detector_buffer& buffer,
                                       callable_t&& callable) {
  return detail::visit_erased_detector(
      buffer,
      [&callable]<concepts::detector detector_t>(
          const detector_buffer_t<detector_t>& concrete) -> decltype(auto) {
        return std::forward<callable_t>(callable)
            .template operator()<detector_t>(get_data(concrete));
      },
      static_cast<detector_list_t*>(nullptr));
}

/// Copy a supported host detector to a buffer using the supplied allocator and
/// copy implementation. Synchronization follows detray::get_buffer's contract.
template <typename detector_list_t>
detector_buffer buffer_from_host_detector(const host_detector& detector,
                                          vecmem::memory_resource& mr,
                                          vecmem::copy& copy) {
  return host_detector_visitor<detector_list_t>(
      detector,
      [&mr, &copy]<concepts::detector detector_t>(const detector_t& concrete) {
        detector_buffer result;
        result.set<detector_t>(get_buffer(concrete, mr, copy));
        return result;
      });
}

}  // namespace detray

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "detray/core/detail/erased_detector_visitor.hpp"
#include "detray/core/detector.hpp"
#include "detray/utils/move_only_any.hpp"

namespace detray {

/// Owning, move-only host detector with its concrete metadata type erased.
/// The memory resource used to construct the detector must outlive this object.
/// No registry dependency: any concrete detector can be stored. Each operation
/// explicitly supplies the list of types for which it has compiled support.
class host_detector {
 public:
  template <typename T>
  void set(detector_host_t<T>&& object) {
    m_obj.set<detector_host_t<T>>(std::move(object));
  }
  template <typename T>
  bool is() const noexcept {
    return m_obj.is<detector_host_t<T>>();
  }
  template <typename T>
  const detector_host_t<T>& as() const {
    return m_obj.as<detector_host_t<T>>();
  }
  const std::type_info& type() const { return m_obj.type(); }
  bool has_value() const noexcept { return m_obj.has_value(); }
  void reset() noexcept { m_obj.reset(); }

 private:
  move_only_any m_obj;
};

/// Dispatch once on the host to a concrete detector. All callable
/// instantiations must return the same type. References are preserved; they
/// borrow the detector.
template <typename detector_list_t, typename callable_t>
decltype(auto) host_detector_visitor(const host_detector& detector,
                                     callable_t&& callable) {
  return detail::visit_erased_detector(detector,
                                       std::forward<callable_t>(callable),
                                       static_cast<detector_list_t*>(nullptr));
}

}  // namespace detray

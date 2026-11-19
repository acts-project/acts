// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include <pybind11/pybind11.h>
#include <vecmem/memory/memory_resource.hpp>

#include "detray/core/host_detector.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/detectors/metadata.hpp"

namespace detray::python {
using algebra = detray::array<DETRAY_CUSTOM_SCALARTYPE>;
using detector_types = registry::detector_types<algebra>;

template <typename callable_t>
decltype(auto) visit(const host_detector& detector, callable_t&& callable) {
  return host_detector_visitor<detector_types>(
      detector, std::forward<callable_t>(callable));
}

/// Retain the Python allocator wrapper with the canonical detector wrapper.
inline pybind11::object wrap(
    host_detector&& detector,
    std::shared_ptr<vecmem::memory_resource> resource) {
  auto object =
      pybind11::cast(std::make_shared<host_detector>(std::move(detector)));
  pybind11::detail::keep_alive_impl(object,
                                    pybind11::cast(std::move(resource)));
  return object;
}

template <concepts::detector detector_t>
pybind11::object wrap(detector_t&& concrete,
                      std::shared_ptr<vecmem::memory_resource> resource) {
  host_detector detector;
  detector.set<detector_t>(std::move(concrete));
  return wrap(std::move(detector), std::move(resource));
}
}  // namespace detray::python

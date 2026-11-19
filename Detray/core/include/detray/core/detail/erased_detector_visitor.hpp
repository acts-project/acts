// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include "detray/core/concepts.hpp"

namespace detray::detail {

template <typename erased_t, typename callable_t>
void visit_erased_detector(const erased_t&, callable_t&&, std::tuple<>*) {
  throw std::invalid_argument("No detector types enabled for this operation");
}

template <typename erased_t, typename callable_t, concepts::detector detector_t,
          concepts::detector... detector_ts>
decltype(auto) visit_erased_detector(const erased_t& object,
                                     callable_t&& callable,
                                     std::tuple<detector_t, detector_ts...>*) {
  if (!object.has_value()) {
    throw std::invalid_argument("Cannot visit an empty detector");
  }
  if (object.template is<detector_t>()) {
    return std::forward<callable_t>(callable).template operator()<detector_t>(
        object.template as<detector_t>());
  }
  if constexpr (sizeof...(detector_ts) > 0) {
    return visit_erased_detector(
        object, std::forward<callable_t>(callable),
        static_cast<std::tuple<detector_ts...>*>(nullptr));
  } else {
    throw std::invalid_argument(std::string("Detector type '") +
                                object.type().name() +
                                "' is not enabled for this operation");
  }
}

}  // namespace detray::detail

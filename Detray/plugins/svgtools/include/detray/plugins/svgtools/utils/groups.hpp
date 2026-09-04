// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <string>
#include <vector>

namespace detray::svgtools::utils {

template <typename iterator_t>
inline auto group(const std::string& identification,
                  const iterator_t& iterator) {
  actsvg::svg::object ret;
  ret._tag = "g";
  ret._id = identification;
  for (const auto& item : iterator) {
    ret.add_object(item);
  }
  return ret;
}

inline auto group(const std::string& identification) {
  return group(identification, std::vector<actsvg::svg::object>{});
}
}  // namespace detray::svgtools::utils

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/io/utils/file_handle.hpp"

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <ios>
#include <string>
#include <vector>

namespace detray::svgtools {

/// @brief Writes a collection of svgs objects to a single file.
template <typename container_t>
inline void write_svg(const std::string& path, const container_t& svgs,
                      bool replace = true) {
  actsvg::svg::file file;
  for (const actsvg::svg::object& obj : svgs) {
    file.add_object(obj);
  }
  std::ios_base::openmode io_mode =
      replace ? std::ios::out | std::ios::trunc : std::ios::out;

  detray::io::file_handle stream{path, ".svg", io_mode};
  *stream << file;
}

/// @brief Writes an svg object to a file.
/// @note To avoid conflict, the ids of the svg objects must be unique.
inline void write_svg(const std::string& path, const actsvg::svg::object& svg,
                      bool replace = true) {
  write_svg(path, std::array{svg}, replace);
}

/// @brief Writes an svg objects to a file.
/// @note To avoid conflict, the ids of the svg objects must be unique.
inline void write_svg(const std::string& path,
                      const std::initializer_list<actsvg::svg::object>& svgs,
                      bool replace = true) {
  std::vector<actsvg::svg::object> arg = svgs;
  write_svg(path, arg, replace);
}

}  // namespace detray::svgtools

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/concepts.hpp"
#include "detray/io/backend/geometry_reader.hpp"
#include "detray/io/backend/homogeneous_material_reader.hpp"
#include "detray/io/backend/material_map_reader.hpp"
#include "detray/io/backend/surface_grid_reader.hpp"
#include "detray/io/frontend/detail/detector_components_reader.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/json/json_converter.hpp"
#include "detray/io/utils/io_metadata.hpp"
#include "detray/utils/logging.hpp"

// System include(s)
#include <filesystem>
#include <ios>
#include <stdexcept>
#include <string>
#include <vector>

namespace detray::io::detail {

/// @brief Function that reads the common header part of a file in json format
inline common_header_payload deserialize_json_header(
    const std::string& file_name) {
  // Read json file
  io::file_handle file{file_name, std::ios_base::in | std::ios_base::binary};
  nlohmann::json in_json;
  *file >> in_json;

  // Reads the header from file
  header_payload<> h = in_json["header"];

  // Need only the common part here
  const common_header_payload& header = h.common;

  if (header.tag < io::detail::minimal_io_version) {
    DETRAY_WARN_HOST("File was generated with a different detray version");
  }

  return header;
}

/// From the list of files that are given @param files, infer the readers that
/// are needed by peeking into the file headers
///
/// @tparam CAP surface grid bin capacity (@TODO make runtime)
/// @tparam DIM dimension of the surface grids, usually 2D
/// @tparam detector_t type of the detector instance: Must match the data that
///                    is read from file!
template <std::size_t CAP, std::size_t DIM, class detector_t>
inline void add_json_readers(
    io::detail::detector_components_reader<detector_t>& reader,
    const std::vector<std::string>& files) noexcept(false) {
  for (const std::filesystem::path file_name : files) {
    if (file_name.empty()) {
      DETRAY_WARN_HOST("Empty file name. Component will not be built");
      continue;
    }

    // Only add readers for json files
    if (file_name.extension() != ".json") {
      continue;
    }

    // Peek at the header to determine the kind of reader that is needed
    auto header = deserialize_json_header(file_name);

    if (header.tag == "geometry") {
      reader.set_detector_name(header.detector);

      using json_geometry_reader = json_converter<detector_t, geometry_reader>;

      reader.template add_geometry_reader<json_geometry_reader>(
          static_cast<std::string>(file_name));

    } else if (header.tag == "homogeneous_material") {
      if constexpr (detray::concepts::has_homogeneous_material<detector_t>) {
        using json_hom_material_reader =
            json_converter<detector_t, homogeneous_material_reader>;

        reader.template add<json_hom_material_reader>(file_name);
      } else {
        print_type_warning<detector_t>(header.tag);
      }
    } else if (header.tag == "material_maps") {
      if constexpr (detray::concepts::has_material_maps<detector_t>) {
        using json_material_map_reader = json_converter<
            detector_t,
            material_map_reader<std::integral_constant<std::size_t, DIM>>>;

        reader.template add<json_material_map_reader>(file_name);
      } else {
        print_type_warning<detector_t>(header.tag);
      }
    } else if (header.tag == "surface_grids") {
      if constexpr (detray::concepts::has_surface_grids<detector_t>) {
        using json_surface_grid_reader = json_converter<
            detector_t,
            surface_grid_reader<typename detector_t::surface_type,
                                std::integral_constant<std::size_t, CAP>,
                                std::integral_constant<std::size_t, DIM>>>;

        reader.template add<json_surface_grid_reader>(file_name);
      } else {
        print_type_warning<detector_t>(header.tag);
      }
    } else {
      throw std::invalid_argument("Unsupported file tag '" + header.tag +
                                  "' in input file: " + file_name.string());
    }
  }
}

}  // namespace detray::io::detail

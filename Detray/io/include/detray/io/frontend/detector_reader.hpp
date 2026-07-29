// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/detector_builder.hpp"
#include "detray/io/frontend/detail/detector_components_reader.hpp"
#include "detray/io/frontend/detector_reader_config.hpp"
#include "detray/io/json/json_converter.hpp"
#include "detray/utils/consistency_checker.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/print_detector.hpp"

// System include(s)
#include <filesystem>
#include <ios>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace detray::io {

namespace detail {

/// @brief Register file converters in the detector components reader
///
/// @tparam input_converter_t the type of converter to be registered
/// @tparam detector_t the type of detector under construction
///
/// @param detector_components class that holds all required readers and converters
/// @param file_names list of input files that get associated to the converter by file extension matching
template <typename input_converter_t, typename detector_t, std::size_t CAP,
          std::size_t DIM>
void add_input_file_converters(detail::detector_components_reader<
                                   detector_t, CAP, DIM>& detector_components,
                               const std::vector<std::string>& file_names) {
  for (const std::filesystem::path file_name : file_names) {
    if (file_name.empty()) {
      DETRAY_DEBUG_HOST("Empty file name. Component will not be built");
      continue;
    }

    // Only add readers for json files
    std::string file_ext{file_name.extension()};
    if (file_ext == input_converter_t{}.file_extension()) {
      detector_components.template add_converter<input_converter_t>(file_name);
    }
  }
}

}  // namespace detail

/// @brief Reader function for detray detectors: payload + optional input files
///
/// @tparam detector_t the type of detector to be built
/// @tparam CAP surface grid bin capacity. If CAP is 0, the grid reader builds
///             a grid type with dynamic bin capacity
/// @tparam DIM dimension of the surface grids, usually 2D
/// @tparam volume_builder_t the type of base volume builder to be used
/// @tparam input_converter_ts converter types for the input file formats
///
/// @param resc the memory resource to be used for the detector container allocs
/// @param cfg the detector reader configuration
/// @param payload [pre-filled] detector intermediate data representation
/// @param file_names optional list of additional input file names.
///
/// @returns a complete detector object + a map that contains the volume names
template <class detector_t, std::size_t CAP = 0u, std::size_t DIM = 2u,
          template <typename> class volume_builder_t = volume_builder,
          typename... input_converter_ts>
auto read_detector(vecmem::memory_resource& resc,
                   const detector_reader_config& cfg,
                   detector_payload& payload) noexcept(false) {
  // Hold all required readers (one for every component) and optional converters
  detail::detector_components_reader<detector_t, CAP, DIM> detector_components;
  const std::vector<std::string>& file_names = cfg.files();

  if constexpr (sizeof...(input_converter_ts) > 0) {
    // Check if data from any input files needs to be converted
    if (!file_names.empty()) {
      // Register the readers for the files in json format
      (add_input_file_converters<input_converter_ts>(detector_components,
                                                     file_names),
       ...);

      // Make sure that all files will be read
      if (detector_components.converter_map().size() != file_names.size()) {
        std::stringstream err_str{};
        for (const auto& conv_pair : detector_components.converter_map()) {
          err_str << "-> " << conv_pair.first << std::endl;
        }
        DETRAY_ERROR_HOST("Not all files were registered to a converter. "
                          << "Successfully registered files:\n"
                          << err_str.str());
      }
    }
  } else {
    assert(file_names.empty());
  }

  // Collect all input data about the detector (from converters and external
  // payload)
  detector_builder<typename detector_t::metadata, volume_builder_t> det_builder;
  detector_components.read(det_builder, payload);

  // Build and return the detector
  auto det = det_builder.build(resc, payload.names);

  if (cfg.do_check()) {
    // This will throw an exception in case of inconsistencies
    detray::detail::check_consistency(det, cfg.verbose_check(), payload.names);
  }

  DETRAY_DEBUG_HOST("\n" << detray::utils::print_detector(det, payload.names));

  return std::make_pair(std::move(det), std::move(payload.names));
}

/// @brief Read the detector completely from input files.
///
/// @tparam detector_t the type of detector to be built
/// @tparam CAP surface grid bin capacity. If CAP is 0, the grid reader builds
///             a grid type with dynamic bin capacity
/// @tparam DIM dimension of the surface grids, usually 2D
/// @tparam volume_builder_t the type of base volume builder to be used
/// @tparam input_converter_ts converter types for the input file formats
///
/// @param resc the memory resource to be used for the detector container allocs
/// @param cfg the detector reader configuration
/// @param file_names list of input file names. Needs at least a geometry file!
///
/// @returns a complete detector object + a map that contains the volume names
template <class detector_t, std::size_t CAP = 0u, std::size_t DIM = 2u,
          template <typename> class volume_builder_t = volume_builder,
          typename... input_converter_ts>
auto read_detector(vecmem::memory_resource& resc,
                   const detector_reader_config& cfg) noexcept(false) {
  // static_assert(sizeof...(input_converter_ts) > 0);

  // Empty payload: All data comes from input files
  detector_payload payload{};

  return read_detector<detector_t, CAP, DIM, volume_builder_t,
                       json_input_converter>(resc, cfg, payload);
}

}  // namespace detray::io

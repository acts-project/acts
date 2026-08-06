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
#include "detray/io/frontend/payloads.hpp"
#include "detray/utils/consistency_checker.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/print_detector.hpp"

// System include(s)
#include <filesystem>
#include <initializer_list>
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
/// @param detector_components class that holds all required readers
///                            and converters
/// @param file_names list of input files that get associated to the converter
///                   by file extension matching
template <typename input_converter_t>
void add_input_file_converters(
    detray::io::detail::detector_components_converter& payload_converter,
    const std::vector<std::string>& file_names) {
  // Register each file according to the file extension
  for (const std::filesystem::path file_name : file_names) {
    if (file_name.empty()) {
      DETRAY_DEBUG_HOST("Empty file name. Component will not be built");
      continue;
    }

    // Only add readers for json files
    std::string file_ext{file_name.extension()};
    if (file_ext == input_converter_t{}.file_extension()) {
      payload_converter.template add_converter<input_converter_t>(file_name);
    }
  }
}

}  // namespace detail

/// @brief Convert input file data into the detector payload
///
/// @param [out] payload the full detector payload, to which the converted data
///                      is added
/// @param [in] file_names list of input file names whose data is converted
///                        into payloads. They are automatically assigned
///                        to a converter via file extensions.
/// @{
template <typename... input_converter_ts>
void convert_to_payload(detray::io::detector_payload& payload,
                        const std::vector<std::string>& file_names) {
  static_assert(sizeof...(input_converter_ts) > 0);

  // Hold all required readers (one for every component) and optional converters
  detray::io::detail::detector_components_converter payload_converter{};

  // Check if data from any input files needs to be converted
  if (!file_names.empty()) {
    // Register the readers for the files in json format
    (detray::io::detail::add_input_file_converters<input_converter_ts>(
         payload_converter, file_names),
     ...);

    // Make sure that all files will be read
    if (payload_converter.converter_map().size() != file_names.size()) {
      std::stringstream err_str{};
      for (const auto& conv_pair : payload_converter.converter_map()) {
        err_str << "-> " << conv_pair.first << std::endl;
      }
      DETRAY_ERROR_HOST(
          "Not all files were registered to a converter: Please check that the "
          "file extensions match the converter types!"
          << "Successfully registered files:\n"
          << err_str.str());
    }
  }

  payload_converter.convert(payload);
}

template <typename... input_converter_ts>
void convert_to_payload(detray::io::detector_payload& payload,
                        std::initializer_list<std::string_view> file_names) {
  detray::io::convert_to_payload<input_converter_ts...>(
      payload,
      std::vector<std::string>(std::begin(file_names), std::end(file_names)));
}
/// @}

/// @brief Reader function for detray detectors: payload + optional input files
///
/// @tparam detector_t the type of detector to be built
/// @tparam GCAP surface grid bin capacity. If CAP is 0, the grid reader builds
///              a grid type with dynamic bin capacity
/// @tparam GDIM dimension of the surface grids, usually 2D
/// @tparam MDIM dimension of the material grids, usually 2D
/// @tparam volume_builder_t the type of base volume builder to be used
/// @tparam input_converter_ts converter types for the input file formats
///
/// @param resc the memory resource to be used for the detector container allocs
/// @param cfg the detector reader configuration
/// @param[in, out] payload [pre-filled] detector intermediate data
///                         representation
///
/// @returns a complete detector object + a map that contains the volume names
template <class detector_t, std::size_t GCAP = 0u, std::size_t GDIM = 2u,
          std::size_t MDIM = 2u,
          template <typename> class volume_builder_t = volume_builder,
          typename... input_converter_ts>
auto read_detector(vecmem::memory_resource& resc,
                   const detray::io::detector_reader_config& cfg,
                   detray::io::detector_payload& payload) noexcept(false) {
  // Convert the input file data and add it to the detector payload object
  if constexpr (sizeof...(input_converter_ts) > 0) {
    detray::io::convert_to_payload<input_converter_ts...>(payload, cfg.files());
  } else {
    // No converters required
    assert(cfg.files().empty());
  }

  // Collect all input data about the detector in the detector builder
  // (from converters and external payload)
  detray::io::detail::detector_components_reader<detector_t, GCAP, GDIM, MDIM>
      detector_components;
  detray::detector_builder<typename detector_t::metadata, volume_builder_t>
      det_builder;

  // Transcribe the data into the detector builder
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
/// @tparam GCAP surface grid bin capacity. If CAP is 0, the grid reader builds
///              a grid type with dynamic bin capacity
/// @tparam GDIM dimension of the surface grids, usually 2D
/// @tparam MDIM dimension of the material grids, usually 2D
/// @tparam volume_builder_t the type of base volume builder to be used
/// @tparam input_converter_ts converter types for the input file formats
///
/// @param resc the memory resource to be used for the detector container allocs
/// @param cfg the detector reader configuration
/// @param file_names list of input file names. Needs at least a geometry file!
///
/// @returns a complete detector object + a map that contains the volume names
template <class detector_t, std::size_t GCAP = 0u, std::size_t GDIM = 2u,
          std::size_t MDIM = 2u,
          template <typename> class volume_builder_t = volume_builder,
          typename... input_converter_ts>
auto read_detector(
    vecmem::memory_resource& resc,
    const detray::io::detector_reader_config& cfg) noexcept(false) {
  // Empty payload: All data comes from input files
  detray::io::detector_payload payload{};

  return detray::io::read_detector<detector_t, GCAP, GDIM, MDIM,
                                   volume_builder_t, input_converter_ts...>(
      resc, cfg, payload);
}

}  // namespace detray::io

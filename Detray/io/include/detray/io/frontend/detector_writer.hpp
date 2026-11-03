// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/io/backend/detail/basic_converter.hpp"
#include "detray/io/frontend/detail/detector_components_writer.hpp"
#include "detray/io/frontend/detector_writer_config.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/utils/create_path.hpp"
#include "detray/utils/logging.hpp"

// System include(s)
#include <filesystem>
#include <ios>
#include <stdexcept>
#include <string>

namespace detray::io {

namespace detail {

/// @brief Register file converters in the detector components output converter
///
/// @tparam output_converter_t the type of converter to be registered
///
/// @param payload_converter class that holds all required output converters
/// @param cfg the detector writer configuration. The converter is only
///            registered, if it matches the configured output format
/// @param file_path the directory to which the converter writes its files
/// @param mode the mode in which the converter opens the files
template <typename output_converter_t>
void add_output_file_converters(
    detray::io::detail::detector_components_output_converter& payload_converter,
    const detray::io::detector_writer_config& cfg,
    const std::filesystem::path& file_path,
    const std::ios_base::openmode mode) {
  // Only add converters for the configured output format
  if (output_converter_t{}.format() == cfg.format()) {
    payload_converter.template add_converter<output_converter_t>(file_path,
                                                                 mode);
  }
}

}  // namespace detail

/// @brief Convert the detector payload into output files
///
/// @tparam output_converter_ts converter types for the output file formats
///
/// @param [in, out] payload the full detector payload, which is converted.
///                          Its header is generated from the payload data
/// @param [in] cfg the detector writer configuration
template <typename... output_converter_ts>
void convert_from_payload(detray::io::detector_payload& payload,
                          const detray::io::detector_writer_config& cfg) {
  static_assert(sizeof...(output_converter_ts) > 0);

  // Check that the required geometry data is there
  if (payload.geometry.volumes.empty()) {
    const std::string err_msg{"No geometry data found in detector payload"};
    DETRAY_FATAL_HOST(err_msg);
    throw std::runtime_error(err_msg);
  }

  // In case the conversion is started from an external payload:
  // Fill the sub-headers (no further detector data available at this level)
  detail::fill_sub_header_payloads(payload, cfg);
  payload.header.source = cfg.source();

  // How to open the file
  const std::ios_base::openmode out_mode{std::ios_base::out |
                                         std::ios_base::binary};
  const auto mode =
      cfg.replace_files() ? (out_mode | std::ios_base::trunc) : out_mode;

  const auto file_path = detray::io::create_path(cfg.path());

  if (cfg.compactify_json()) {
    DETRAY_ERROR_HOST("Compactifying json files is not yet implemented");
  }

  // Hold all required output converters
  detray::io::detail::detector_components_output_converter payload_converter{};

  // Register the converters that match the output format
  (detray::io::detail::add_output_file_converters<output_converter_ts>(
       payload_converter, cfg, file_path, mode),
   ...);

  // Make sure that the detector data will be written
  if (payload_converter.converters().empty()) {
    const std::string err_msg{
        "No converter registered for the configured output format: Please "
        "check that the output format matches the converter types!"};
    DETRAY_FATAL_HOST(err_msg);
    throw std::invalid_argument(err_msg);
  }

  payload_converter.convert(payload);
}

/// @brief Writer function for detray detectors: some components are written
/// from externally provided payload
///
/// Based on both the given config/file format and the detector type,
/// the correct writers and converters are being assembled and called.
///
/// @tparam detector_t the type of detector to be written
/// @tparam output_converter_ts converter types for the output file formats
///
/// @param det the detector to be written
/// @param names the volume names of the detector
/// @param cfg the detector writer configuration
/// @param[out] payload detector intermediate data representation
template <class detector_t, typename... output_converter_ts>
void write_detector(const detector_t& det,
                    const typename detector_t::name_map& names,
                    detray::io::detector_writer_config& cfg,
                    detray::io::detector_payload& payload) noexcept(false) {
  DETRAY_INFO_HOST("Writing detray detector...");

  if (det.material_store().all_empty() &&
      !payload.homogeneous_material.has_value() &&
      !payload.material_maps.has_value()) {
    cfg.write_material(false);
  }

  // Collect all data about the detector in the detector payload
  detray::io::detail::detector_components_writer<detector_t>
      detector_components;
  detector_components.write(det, names, cfg, payload);

  // Convert the detector payload to the output format
  if constexpr (sizeof...(output_converter_ts) > 0) {
    detray::io::convert_from_payload<output_converter_ts...>(payload, cfg);
  } else {
    DETRAY_INFO_HOST(
        "No output converters given to 'detector_writer'. Only the payload "
        "will be filled...");
  }

  DETRAY_INFO_HOST("Done writing: " << det.name(names)
                                    << " (Source: " << cfg.source() << ")");
}

/// @brief Write the output file completely from detector
///
/// @tparam detector_t the type of detector to be written
/// @tparam output_converter_ts converter types for the output file formats
///
/// @param det the detector to be written
/// @param names the volume names of the detector
/// @param cfg the detector writer configuration
template <class detector_t, typename... output_converter_ts>
void write_detector(const detector_t& det,
                    const typename detector_t::name_map& names,
                    detray::io::detector_writer_config& cfg) noexcept(false) {
  // No converters present: Output files cannot be written
  if constexpr (sizeof...(output_converter_ts) == 0) {
    const std::string err_msg{
        "No converters compiled into the 'write_detector' function: Cannot "
        "write output files!\nUse e.g. the 'write_detector_json' function for "
        "json files or add custom output converter types to 'write_detector'."};
    DETRAY_FATAL_HOST(err_msg);
    throw std::invalid_argument(err_msg);
  }

  // Empty payload to be filled: All data comes from the detector
  detray::io::detector_payload payload{};

  detray::io::write_detector<detector_t, output_converter_ts...>(det, names,
                                                                 cfg, payload);
}

}  // namespace detray::io

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
#include "detray/io/backend/geometry_writer.hpp"
#include "detray/io/backend/homogeneous_material_writer.hpp"
#include "detray/io/backend/material_map_writer.hpp"
#include "detray/io/backend/surface_grid_writer.hpp"
#include "detray/io/frontend/detector_writer_config.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/frontend/writer_interface.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/type_list.hpp"

// System include(s)
#include <any>
#include <filesystem>
#include <ios>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

namespace detray::io::detail {

/// @brief Convert a detector payload to output data
///
/// The class aggregates a number of different converters depending on the
/// output sink and calls them once the detector data should be written
/// (e.g to output files).
class detector_components_output_converter {
  using converter_ptr_t = std::unique_ptr<output_converter_interface>;
  using file_converter_ptr_t = std::unique_ptr<output_file_converter_interface>;

 public:
  /// Default constructor
  detector_components_output_converter() = default;

  /// Create a new converter of type @tparam converter_t
  ///
  /// @param file_path the directory in which the converter writes its files
  /// @param mode the mode in which the converter opens the files
  template <class converter_t>
    requires std::is_base_of_v<output_file_converter_interface, converter_t>
  void add_converter(const std::filesystem::path& file_path,
                     const std::ios_base::openmode mode) {
    auto conv_ptr = std::make_unique<converter_t>();
    conv_ptr->open_mode(mode);
    add_converter(std::move(conv_ptr), file_path);
  }

  /// Attach an existing converter via @param conv_ptr
  void add_converter(file_converter_ptr_t&& conv_ptr,
                     const std::filesystem::path& file_path) {
    m_converters.emplace_back(file_path, std::move(conv_ptr));
  }

  /// @returns access to the converters - const
  const auto& converters() const { return m_converters; }

  /// Run all registered converters for their corresponding output sink
  ///
  /// @param[in] payload detector data rep. that is converted
  void convert(const detray::io::detector_payload& payload) {
    // Nothing left to do
    if (m_converters.empty()) {
      return;
    }

    // Convert the detector payload into the output data
    for (const auto& [file_path, converter] : m_converters) {
      converter->from_payload(payload, std::any{file_path});
    }

    DETRAY_INFO_HOST(" -> Done converting output data");
  }

 private:
  std::vector<std::pair<std::filesystem::path, converter_ptr_t>> m_converters;
};

/// @brief Writes the detector data into the detector payload
///
/// Automatically detects the kind of writers that are needed from the
/// detector type and the writer configuration
template <class detector_t>
class detector_components_writer final {
  using writer_ptr_t = std::unique_ptr<writer_interface<detector_t>>;

 public:
  /// Default constructor
  detector_components_writer() = default;

  /// Writes the full detector data of @param det into the @param payload by
  /// calling the writers, while using the name map @param names for the
  /// volume names.
  ///
  /// @param[in] det the detector that is converted
  /// @param[in] names the volume names of the detector
  /// @param[in] cfg the detector writer configuration: selects the components
  /// @param[out] payload detector data rep. filled by the writers
  void write(const detector_t& det, const typename detector_t::name_map& names,
             const detray::io::detector_writer_config& cfg,
             detray::io::detector_payload& payload) {
    // Geometry writer: Always required
    std::vector<writer_ptr_t> writers;
    writers.push_back(
        std::make_unique<detray::io::geometry_writer<detector_t>>());

    // Add any optional writers
    if (cfg.write_material()) {
      if constexpr (detray::concepts::has_homogeneous_material<detector_t>) {
        writers.push_back(
            std::make_unique<
                detray::io::homogeneous_material_writer<detector_t>>());
      }
      if constexpr (detray::concepts::has_material_maps<detector_t>) {
        writers.push_back(
            std::make_unique<detray::io::material_map_writer<detector_t>>());
      }
    }
    if constexpr (detray::concepts::has_surface_grids<detector_t>) {
      if (cfg.write_grids()) {
        writers.push_back(
            std::make_unique<detray::io::surface_grid_writer<detector_t>>());
      }
    }

    // Set the detector name and metadata type
    payload.detector_name = det.name(names);
    payload.header.metadata =
        detray::types::get_name<typename detector_t::metadata>();

    // Write the detector components into the payload
    for (const auto& writer : writers) {
      writer->to_payload(det, names, payload);
    }
    DETRAY_INFO_HOST(" -> Done writing payloads");
  }
};

}  // namespace detray::io::detail

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
#include "detray/builders/volume_builder.hpp"
#include "detray/io/backend/geometry_reader.hpp"
#include "detray/io/backend/homogeneous_material_reader.hpp"
#include "detray/io/backend/material_map_reader.hpp"
#include "detray/io/backend/surface_grid_reader.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/frontend/reader_interface.hpp"

// System include(s)
#include <algorithm>
#include <cassert>
#include <map>
#include <memory>
#include <string>
#include <type_traits>

namespace detray::io::detail {

/// @brief Convert input data to a detector payload
///
/// The class aggregates a number of different converters depending on the
/// input source and calls them once the detector data should be read in
/// (e.g from input files).
class detector_components_converter {
  using converter_ptr_t = std::unique_ptr<input_converter_interface>;
  using file_converter_ptr_t = std::unique_ptr<input_file_converter_interface>;

 public:
  /// Default constructor
  detector_components_converter() = default;

  /// Create a new converter of type @tparam converter_t
  template <class converter_t>
    requires std::is_base_of_v<input_file_converter_interface, converter_t>
  void add_converter(const std::string& file_name) {
    add_converter(std::make_unique<converter_t>(), file_name);
  }

  /// Attach an existing converter via @param conv_ptr
  void add_converter(file_converter_ptr_t&& conv_ptr,
                     const std::string& file_name) {
    m_converters[file_name] = std::move(conv_ptr);
  }

  /// @returns access to the converters map - const
  const auto& converter_map() const { return m_converters; }

  /// Run all registered converters for their corresponding input data
  ///
  /// @param[out] payload detector data rep. filled by converters
  void convert(detray::io::detector_payload& payload) {
    // Nothing left to do
    if (m_converters.empty()) {
      return;
    }

    // Convert the requested detector components into their payloads and add
    // them to the detector_payload
    for (const auto& [file_name, converter] : m_converters) {
      converter->add_to_payload(file_name, payload);
    }

    DETRAY_INFO_HOST(" -> Done converting input data");
  }

 private:
  std::unordered_map<std::string, converter_ptr_t> m_converters;
};

/// @brief Reads the detector payload into the detector builder
///
/// Automatically detects the kind of readers that are needed to process
/// the payload object at runtime
template <class detector_t, std::size_t GCAP = 0u, std::size_t GDIM = 2u,
          std::size_t MDIM = 2u>
class detector_components_reader final {
  using reader_ptr_t = std::unique_ptr<reader_interface<detector_t>>;

 public:
  /// Default constructor
  detector_components_reader() = default;

  /// Reads the full detector into @param det by calling the readers, while
  /// using the name map @param volume_names for to write the volume names.
  ///
  /// @param[out] det_builder complete detector data filled by readers
  ///                         from payload
  /// @param[in] payload externally provided detector data, filled into
  ///                    the detector builder by readers held in this class
  void read(detray::detector_builder<typename detector_t::metadata,
                                     volume_builder>& det_builder,
            detray::io::detector_payload& payload) {
    // Check that the required geometry data is there
    if (payload.geometry.volumes.empty()) {
      std::string err_msg{"No geometry data found in input"};
      DETRAY_FATAL_HOST(err_msg);
      throw std::runtime_error(err_msg);
    }

    // Geometry reader: Always required
    auto geo_reader =
        std::make_unique<detray::io::geometry_reader<detector_t>>();

    // Add any optional readers
    std::vector<reader_ptr_t> readers;
    if constexpr (detray::concepts::has_surface_grids<detector_t>) {
      if (payload.surface_grids.has_value()) {
        using GRID_CAP = std::integral_constant<std::size_t, GCAP>;
        using GRID_DIM = std::integral_constant<std::size_t, GDIM>;

        readers.push_back(
            std::make_unique<
                surface_grid_reader<detector_t, GRID_CAP, GRID_DIM>>());
      }
    }
    if constexpr (detray::concepts::has_homogeneous_material<detector_t>) {
      if (payload.homogeneous_material.has_value()) {
        readers.push_back(
            std::make_unique<
                detray::io::homogeneous_material_reader<detector_t>>());
      }
    }
    if constexpr (detray::concepts::has_material_maps<detector_t>) {
      if (payload.material_maps.has_value()) {
        using MAP_DIM = std::integral_constant<std::size_t, MDIM>;
        readers.push_back(
            std::make_unique<
                detray::io::material_map_reader<detector_t, MAP_DIM>>());
      }
    }

    // Set the detector name (may have been read from input file)
    det_builder.set_name(payload.detector_name);

    // Fill the geometry data first
    geo_reader->from_payload(det_builder, payload);

    // Read the other component payloads into the detector builder
    for (const auto& reader : readers) {
      reader->from_payload(det_builder, payload);
    }
    DETRAY_INFO_HOST(" -> Done reading payloads");
  }
};

}  // namespace detray::io::detail

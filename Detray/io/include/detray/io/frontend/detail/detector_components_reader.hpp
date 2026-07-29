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
#include "detray/io/frontend/reader_interface.hpp"

// System include(s)
#include <algorithm>
#include <cassert>
#include <map>
#include <memory>
#include <string>
#include <type_traits>

namespace detray::io::detail {

/// @brief A reader that automatically detects the kind of backend readers and converters to load in order to
///
/// The class aggregates a number of different readers and calls them once
/// the detector data should be read in from file.
template <class detector_t, std::size_t CAP = 0u, std::size_t DIM = 2u>
class detector_components_reader final {
  using reader_ptr_t = std::unique_ptr<reader_interface<detector_t>>;
  using converter_ptr_t = std::unique_ptr<input_converter_interface>;

 public:
  /// Default constructor
  detector_components_reader() = default;

  /// Create a new converter of type @tparam converter_t
  template <class converter_t>
    requires std::is_base_of_v<input_file_converter_interface, converter_t>
  void add_converter(const std::string& file_name) {
    add_converter(std::make_unique<converter_t>(), file_name);
  }

  /// Attach an existing reader via @param r_ptr to the readers
  void add_converter(converter_ptr_t&& r_ptr, const std::string& file_name) {
    m_converters[file_name] = std::move(r_ptr);
  }

  /// @returns access to the converters map - const
  const auto& converter_map() const { return m_converters; }

  /// @returns access to the readers vector - const
  const auto& readers_map() const { return m_readers; }

  /// Run all registered converters for their corresponding input data
  ///
  /// @param[out] payload detector data rep. filled by converters
  void convert(detector_payload& payload) {
    // Nothing left to do
    if (m_converters.empty()) {
      return;
    }

    // Read the other component payloads into the detector builder
    for (const auto& [file_name, converter] : m_converters) {
      converter->convert_to_payload(file_name, payload);
    }

    DETRAY_INFO_HOST(" -> Done converting input data");
  }

  /// Reads the full detector into @param det by calling the readers, while
  /// using the name map @param volume_names for to write the volume names.
  ///
  /// @param[out] det_builder complete detector data filled by readers from payload
  /// @param[in, out] payload externally provided detector data, additionally filled by converters registered in this class
  void read(detector_builder<typename detector_t::metadata, volume_builder>&
                det_builder,
            detector_payload& payload) {
    // If there are any input converters registered, let them fill the payload
    convert(payload);

    // Once all data is ready, dynamically register the corresponding readers
    if (payload.geometry.volumes.empty()) {
      throw std::runtime_error("No geometry data found in detector reader!");
    }

    m_geo_reader = std::make_unique<geometry_reader<detector_t>>();

    if constexpr (concepts::has_surface_grids<detector_t>) {
      if (payload.surface_grids.has_value()) {
        // Dynamic size of bin content
        using GRID_CAP = std::integral_constant<std::size_t, CAP>;
        // 2D grids
        using GRID_DIM = std::integral_constant<std::size_t, DIM>;

        m_readers.push_back(
            std::make_unique<
                surface_grid_reader<detector_t, GRID_CAP, GRID_DIM>>());
      }
    }
    if constexpr (concepts::has_homogeneous_material<detector_t>) {
      if (payload.homogeneous_material.has_value()) {
        m_readers.push_back(
            std::make_unique<homogeneous_material_reader<detector_t>>());
      }
    }
    if constexpr (concepts::has_material_maps<detector_t>) {
      if (payload.material_maps.has_value()) {
        // 2D material maps
        using MAP_DIM = std::integral_constant<std::size_t, 2u>;
        m_readers.push_back(
            std::make_unique<material_map_reader<detector_t, MAP_DIM>>());
      }
    }

    // Set the detector name (may have been read from input file)
    det_builder.set_name(payload.names.get_detector_name());

    // We have to at least read a geometry
    if (m_geo_reader == nullptr) {
      throw std::runtime_error("No geometry reader registered!");
    }
    // Fill the geometry data first
    m_geo_reader->read_from_payload(det_builder, payload);

    // Read the other component payloads into the detector builder
    for (const auto& reader : m_readers) {
      reader->read_from_payload(det_builder, payload);
    }
    DETRAY_INFO_HOST(" -> Done reading payloads");
  }

 private:
  /// Payload reader that fills the geometry into the detector builder
  reader_ptr_t m_geo_reader{nullptr};
  /// The backend readers registered to fill the detector builder: e.g.
  /// material, grids...
  std::vector<reader_ptr_t> m_readers;
  /// The input converters registered to generate the payloads: e.g. material,
  /// grids...
  std::unordered_map<std::string, converter_ptr_t> m_converters;
};

}  // namespace detray::io::detail

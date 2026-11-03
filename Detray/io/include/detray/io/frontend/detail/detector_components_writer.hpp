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

namespace detray::io {

namespace detail {

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

}  // namespace detail

/// @brief Generate common header information for a given detector
/// @tparam detector_t the detector type
/// @param det The detector for which to generate the header data
/// @param source Description of the source of the detector
/// @param names Map of volume names
/// @returns the detector header payload with data common to all components
template <detray::concepts::detector detector_t>
inline common_header_payload generate_common_header(
    const detector_t& det, const std::string& source,
    const typename detector_t::name_map& names = {}) {
  // Once the payload data is filled, generate and add the header information
  common_header_payload header{};

  header.version = io::detail::get_detray_version();
  header.date = io::detail::get_current_date();
  header.metadata = detray::types::get_name<typename detector_t::metadata>();
  header.source = source;

  if (!names.empty()) {
    header.detector = det.name(names);
  }

  return header;
}

namespace detail {

/// @brief Fill the sub-header part for every detector component
/// @param [in, out] payload The detector payload to which the header data is added
/// @param cfg The detector writer configuration
inline void fill_sub_header_payloads(
    detector_payload& payload,
    const detray::io::detector_writer_config& cfg = {}) {
  // Once the payload data is filled, generate and add the header information
  header_payload& header = payload.header;

  const auto has_data =
      []<typename payload_t>(
          const std::optional<payload_t>& sub_payload) -> bool {
    return sub_payload.has_value() && sub_payload->is_valid();
  };

  // If the payload is empty, do nothing
  if (payload.geometry.is_valid()) {
    header.geometry = to_header_payload(payload.geometry);
  }
  if (cfg.write_material()) {
    if (has_data(payload.homogeneous_material)) {
      header.homogeneous_material =
          to_header_payload(*payload.homogeneous_material);
    }
    if (has_data(payload.material_maps)) {
      header.material_maps = to_header_payload(*payload.material_maps);
    }
  }

  if (cfg.write_grids()) {
    if (has_data(payload.surface_grids)) {
      header.surface_grids = to_header_payload(*payload.surface_grids);
    }
    if (has_data(payload.volume_grids)) {
      header.volume_grids = to_header_payload(*payload.volume_grids);
    }
  }
}

/// @brief Fill the header part of a detector payload (including common data)
/// @tparam detector_t the detector type
/// @param det The detector for which to generate the header data
/// @param [in, out] payload The detector payload to which the header data is added
/// @param cfg The detector writer configuration
/// @param names Map of volume names
template <detray::concepts::detector detector_t>
inline void fill_header_payload(
    const detector_t& det, detector_payload& payload,
    const detray::io::detector_writer_config& cfg,
    const typename detector_t::name_map& names = {}) {
  // Fill the common information
  payload.header = generate_common_header(det, cfg.source(), names);

  // Fill the subheaders from the payload data
  fill_sub_header_payloads(payload, cfg);
}

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
    // Always take the data of the actual detector
    const auto print_warning =
        []([[maybe_unused]] const std::string& payload_type) {
          DETRAY_WARN_HOST("Writer found "
                           << payload_type
                           << " data already present in payload: Detector data "
                              "will not be converted!");
        };
    const auto throw_type_err = [](const std::string& payload_type) {
      std::stringstream stream{};
      stream << "Writer found " << payload_type << " payload, but detector ("
             << detray::types::get_name<typename detector_t::metadata>()
             << ") does not support it: Check the metadata!";
      DETRAY_FATAL_HOST(stream.str());
      throw std::invalid_argument(stream.str());
    };

    // Either take the geo. payload or write the detector geometry to payload
    std::vector<writer_ptr_t> writers;
    if (payload.geometry.volumes.empty()) {
      writers.push_back(
          std::make_unique<detray::io::geometry_writer<detector_t>>());
    } else {
      print_warning("geometry");
    }

    // Add any optional writers
    const auto has_data =
        []<typename payload_t>(
            const std::optional<payload_t>& sub_payload) -> bool {
      return sub_payload.has_value() && sub_payload->is_valid();
    };
    if (cfg.write_material()) {
      const bool contains_hom_material{has_data(payload.homogeneous_material)};
      const bool contains_material_maps{has_data(payload.material_maps)};

      if constexpr (detray::concepts::has_homogeneous_material<detector_t>) {
        if (!contains_hom_material) {
          writers.push_back(
              std::make_unique<
                  detray::io::homogeneous_material_writer<detector_t>>());
        } else {
          print_warning("homogeneous material");
        }
      } else if (contains_hom_material) {
        throw_type_err("homogeneous material");
      }

      if constexpr (detray::concepts::has_material_maps<detector_t>) {
        if (!contains_material_maps) {
          writers.push_back(
              std::make_unique<detray::io::material_map_writer<detector_t>>());
        } else {
          print_warning("material maps");
        }
      } else if (contains_material_maps) {
        throw_type_err("material maps");
      }
    }

    if (cfg.write_grids()) {
      const bool contains_sf_grids{has_data(payload.surface_grids)};

      if constexpr (detray::concepts::has_surface_grids<detector_t>) {
        if (!contains_sf_grids) {
          writers.push_back(
              std::make_unique<detray::io::surface_grid_writer<detector_t>>());
        } else {
          print_warning("surface grids");
        }
      } else if (contains_sf_grids) {
        throw_type_err("surface grids");
      }
    }

    // Write the selected detector components into the payload
    for (const auto& writer : writers) {
      writer->to_payload(det, names, payload);
    }

    // Add the header information to the payload
    fill_header_payload(det, payload, cfg, names);

    DETRAY_INFO_HOST(" -> Done writing payloads");
  }
};

}  // namespace detail
}  // namespace detray::io

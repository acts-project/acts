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
#include "detray/io/backend/detail/basic_converter.hpp"
#include "detray/io/frontend/definitions.hpp"
#include "detray/io/frontend/reader_interface.hpp"
#include "detray/io/frontend/writer_interface.hpp"
#include "detray/io/json/json.hpp"
#include "detray/io/json/json_io.hpp"
#include "detray/io/utils/file_handle.hpp"
#include "detray/io/utils/io_metadata.hpp"

// System include(s)
#include <iomanip>
#include <ios>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace detray::io {

/// @brief Class that reads json data from file and converts it to a detector
/// payload (single component per file)
///
/// This class provides the payload data from the json stream and adds it to the
/// sub-payload of the detector component that is tagged in the file
///
/// @note The class fulfills the @c input_file_converter_interface and
/// can be added to the @c detector_components_converter
class json_input_converter final : public input_file_converter_interface {
 public:
  /// Set json file extension
  json_input_converter() : input_file_converter_interface(".json") {}

  /// Reads json detector data from @param file and adds it to a detector
  /// @param payload object
  void add_from_file(io::file_handle& file,
                     detector_payload& payload) override {
    // Reads the data from file and returns the corresponding io payloads
    nlohmann::json in_json;
    *file >> in_json;

    // Peek at the header to determine the kind of payload is in the file
    const tagged_header_payload header = in_json["header"];
    // If not already done, add common header information to detector payload
    if (payload.header.version.empty()) {
      payload.detector_name = header.detector;
      payload.header = header;
    }

    if (header.tag < io::detail::minimal_io_version) {
      DETRAY_WARN_HOST("File was generated with a different detray version");
    }

    nlohmann::json& json_payload = in_json["data"];
    if (header.tag == "geometry") {
      payload.geometry = json_payload;
    } else if (header.tag == "homogeneous_material") {
      payload.homogeneous_material = json_payload;
    } else if (header.tag == "material_maps") {
      payload.material_maps = json_payload;
    } else if (header.tag == "surface_grids") {
      payload.surface_grids = json_payload;
    } else {
      std::stringstream err_stream{};
      err_stream << "Unsupported file tag '" << header.tag << "' in input file";
      DETRAY_FATAL_HOST(err_stream.str());
      throw std::invalid_argument(err_stream.str());
    }
  }
};

/// @brief Class that converts a detector payload to json and writes it to
/// file (single component per file)
///
/// This class writes every sub-payload of the detector payload that contains
/// data into a separate json file, together with the respective file header.
///
/// @note The class fulfills the @c output_file_converter_interface and
/// can be added to the @c detector_components_output_converter
class json_output_converter final : public output_file_converter_interface {
 public:
  /// Set json file format and extension
  json_output_converter()
      : output_file_converter_interface(io::format::json, ".json") {}

  /// Writes the detector components in @param payload to json files
  ///
  /// @param file_path the directory in which to write the files
  /// @param mode the mode in which the files are opened
  ///
  /// @returns the names of the files that were written
  std::vector<std::string> write_to_file(
      const detector_payload& payload, const std::filesystem::path& file_path,
      const std::ios_base::openmode mode = std::ios::out |
                                           std::ios::binary) override {
    // Assert output stream
    assert(((mode == std::ios_base::out) ||
            (mode == (std::ios_base::out | std::ios_base::binary)) ||
            (mode == (std::ios_base::out | std::ios_base::trunc)) ||
            (mode == (std::ios_base::out | std::ios_base::trunc |
                      std::ios_base::binary))) &&
           "Illegal file mode for json writer");

    const header_payload& header = payload.header;

    // Write every detector component that contains data
    std::vector<std::string> file_names{};
    if (!payload.geometry.volumes.empty()) {
      file_names.push_back(write_component(payload, "geometry", header.geometry,
                                           payload.geometry, file_path, mode));
    }
    if (payload.homogeneous_material.has_value()) {
      file_names.push_back(write_component(
          payload, "homogeneous_material", header.homogeneous_material,
          *payload.homogeneous_material, file_path, mode));
    }
    if (payload.material_maps.has_value()) {
      file_names.push_back(
          write_component(payload, "material_maps", header.material_maps,
                          *payload.material_maps, file_path, mode));
    }
    if (payload.surface_grids.has_value()) {
      file_names.push_back(
          write_component(payload, "surface_grids", header.surface_grids,
                          *payload.surface_grids, file_path, mode));
    }

    return file_names;
  }

 private:
  /// Write a single detector component to a json file
  ///
  /// @param payload the detector payload (header information)
  /// @param tag the tag of the component, which is part of the file name
  /// @param sub_header the component header. It is generated from
  ///                   @param data, if it is missing
  /// @param data the component payload data
  /// @param file_path the directory in which to write the file
  /// @param mode the mode in which the file is opened
  ///
  /// @returns the name of the file that was written
  template <typename sub_header_t, typename data_t>
  std::string write_component(const detector_payload& payload,
                              const std::string_view tag,
                              const std::optional<sub_header_t>& sub_header,
                              const data_t& data,
                              const std::filesystem::path& file_path,
                              const std::ios_base::openmode mode) const {
    // Create a new file
    const std::string file_stem{payload.detector_name + "_" + std::string(tag)};
    const std::string file_ext{file_extension()};
    io::file_handle file{file_path / file_stem, file_ext, mode};

    // Common header information with detector name and content tag
    tagged_header_payload tagged_header{};
    static_cast<common_header_payload&>(tagged_header) = payload.header;
    tagged_header.detector = payload.detector_name;
    tagged_header.tag = tag;

    // Write the header and the detector data into the json stream by using
    // the conversion functions defined in "detray/io/json/json_io.hpp"
    nlohmann::ordered_json out_json;
    out_json["header"] =
        sub_header.has_value()
            ? *sub_header
            : io::detail::basic_converter::to_header_payload(data);
    out_json["header"]["common"] = tagged_header;
    out_json["data"] = data;

    // Write to file
    *file << std::setw(4) << out_json << std::endl;

    return file_stem + file_ext;
  }
};

}  // namespace detray::io

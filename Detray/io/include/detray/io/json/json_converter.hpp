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
#include "detray/io/backend/concepts.hpp"
#include "detray/io/frontend/reader_interface.hpp"
#include "detray/io/frontend/writer_interface.hpp"
#include "detray/io/json/json.hpp"
#include "detray/io/json/json_io.hpp"
#include "detray/io/utils/file_handle.hpp"
#include "detray/io/utils/io_metadata.hpp"

// System include(s)
#include <ios>
#include <iostream>
#include <string>

namespace detray::io {

/// @brief Function that reads the common header part of a file in json format
inline common_header_payload deserialize_json_header(io::file_handle& file) {
  // Read json file
  nlohmann::json in_json;
  std::cout << "Deserializing header" << std::endl;
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

/// @brief Function that reads the common header part of a file in json format
inline common_header_payload deserialize_json_header(
    const std::string& file_name) {
  // Open input file
  io::file_handle file{file_name, std::ios_base::in | std::ios_base::binary};

  return deserialize_json_header(file);
}

/// @brief Class that adds json functionality to backend reader types.
///
/// Assemble the json readers from the backend reader types, which handle the
/// volume builders, and this class, which provides the payload data from the
/// json stream. It also includes the respective @c to_json and @c from_json
/// functions for the payloads ("json_serializers").
///
/// @note The resulting reader types will fulfill @c reader_interface
class json_input_converter final : public input_file_converter_interface {
 public:
  /// Set json file extension
  json_input_converter() : input_file_converter_interface(".json") {}

  /// Writes the geometry to file with a given name
  void from_file(io::file_handle& file, detector_payload& payload) override {
    // Reads the data from file and returns the corresponding io payloads
    nlohmann::json in_json;
    *file >> in_json;

    // Peek at the header to determine the kind of payload is in the file
    header_payload<> h = in_json["header"];

    // Need only the common part here
    const common_header_payload& header = h.common;

    if (header.tag < io::detail::minimal_io_version) {
      DETRAY_WARN_HOST("File was generated with a different detray version");
    }

    payload.names.set_detector_name(header.detector);

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

/// @brief Class that adds converts io payloads to json and back
template <class detector_t, class backend_t>
class json_converter {};

/// @brief Class that adds json functionality to backend writer types.
///
/// Assemble the json writers from the backend writer types, which serialize a
/// detector into the io payloads, and this class, which does the file
/// handling and provides the json stream. It also includes the respective
/// @c to_json and @c from_json functions for the payloads ("json_serializers").
///
/// @note The resulting writer types will fulfill @c writer_interface
template <class detector_t, class backend_t>
  requires concepts::writer_backend<detector_t, backend_t>
class json_converter<detector_t, backend_t> final
    : public writer_interface<detector_t> {
  using io_backend = backend_t;

 public:
  /// File gets created with the json file extension
  json_converter() : writer_interface<detector_t>(".json") {}

  /// Writes the geometry to file with a given name
  std::string write(const detector_t& det,
                    const typename detector_t::name_map& names,
                    const std::ios_base::openmode mode = std::ios::out |
                                                         std::ios::binary,
                    const std::filesystem::path& file_path = {"./"}) override {
    // Assert output stream
    assert(((mode == std::ios_base::out) ||
            (mode == (std::ios_base::out | std::ios_base::binary)) ||
            (mode == (std::ios_base::out | std::ios_base::trunc)) ||
            (mode == (std::ios_base::out | std::ios_base::trunc |
                      std::ios_base::binary))) &&
           "Illegal file mode for json writer");

    std::string det_name = det.name(names);

    // Create a new file
    std::string file_stem{det_name + "_" + std::string(io_backend::tag)};
    io::file_handle file{file_path / file_stem, this->file_extension(), mode};

    // Write some general information
    nlohmann::ordered_json out_json;
    out_json["header"] = io_backend::header_to_payload(det, det_name);

    // Write the detector data into the json stream by using the
    // conversion functions defined in "detray/io/json/json_io.hpp"
    out_json["data"] = io_backend::to_payload(det, names);

    // Write to file
    *file << std::setw(4) << out_json << std::endl;

    return file_stem + this->file_extension();
  }
};

}  // namespace detray::io

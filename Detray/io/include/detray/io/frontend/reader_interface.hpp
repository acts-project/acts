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
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/utils/file_handle.hpp"

// System include(s)
#include <any>
#include <filesystem>
#include <ios>
#include <string>
#include <string_view>

namespace detray::io {

/// @brief Abstract base class for detray detector component readers
///
/// Interface for components readers from payload to @c detector_builder
///
/// @tparam detector_t the detector type under construction
template <class detector_t>
class reader_interface {
 public:
  /// Default destructor
  virtual ~reader_interface() = default;

  /// Returns the string tag that describes the type of reader, which is also
  /// present in data file headers to match the IO data
  virtual std::string_view tag() const = 0;

  /// Reads the respective detector component from file. Since the detector
  /// does not keep the volume names, the name map is also passed and
  /// filled.
  virtual void from_payload(
      detector_builder<typename detector_t::metadata, volume_builder>&,
      const detector_payload&) const = 0;
};

/// @brief Abstract base class for detray detector component converters
class input_converter_interface {
 public:
  /// Default destructor
  virtual ~input_converter_interface() = default;

  /// Reads the respective detector component from an input data source.
  virtual void add_to_payload(const std::any&, detector_payload&) = 0;
};

/// @brief Abstract base class for detray detector component file readers
class input_file_converter_interface : public input_converter_interface {
 public:
  /// Constructor setting a file extension
  explicit input_file_converter_interface(const std::string& ext)
      : m_file_extension{ext} {};

  std::string_view file_extension() const { return m_file_extension; }

  /// Reads the respective detector component from file.
  virtual void add_from_file(io::file_handle&, detector_payload&) = 0;

  /// Reads the respective detector component from file.
  void add_to_payload(const std::any& input_data,
                      detector_payload& payload) override {
    const std::string file_name = get_file_name(input_data);
    io::file_handle file{file_name, std::ios_base::in | std::ios_base::binary};

    add_from_file(file, payload);
  };

 private:
  /// @returns the file name that was passed as type-erased @param input_data
  std::string get_file_name(const std::any& input_data) {
    // Read input data from file
    std::string file_name{""};
    try {
      file_name = std::any_cast<std::string>(input_data);
    } catch (const std::bad_any_cast& e) {
      DETRAY_FATAL_HOST(
          "Unknown input data type in "
          << m_file_extension
          << " file converter! Should be 'std::string' (file name)");
      throw;
    }

    return file_name;
  }

  /// Extension that matches the file format of the respective reader
  std::string m_file_extension;
};

}  // namespace detray::io

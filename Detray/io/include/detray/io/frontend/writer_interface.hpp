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
#include "detray/io/frontend/definitions.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/utils/logging.hpp"

// System include(s)
#include <any>
#include <filesystem>
#include <ios>
#include <string>
#include <string_view>
#include <vector>

namespace detray::io {

/// @brief Abstract base class for detray detector component writers
///
/// Interface for component writers from a detector to a @c detector_payload
///
/// @tparam detector_t the detector type that is converted
template <detray::concepts::detector detector_t>
class writer_interface {
 public:
  /// Default destructor
  virtual ~writer_interface() = default;

  /// @returns the string tag that describes the type of writer, which is also
  /// present in data file headers to match the IO data
  virtual std::string_view tag() const = 0;

  /// Writes the respective detector component to its payload. Since the
  /// detector does not provide the volume names, the name map is also passed.
  /// @param detector detector components to serialize
  /// @param names map used to name detector components
  /// @param payload the payload to be filled with the detector component.
  virtual void to_payload(const detector_t& detector,
                          const typename detector_t::name_map& names,
                          detector_payload& payload) const = 0;
};

/// @brief Abstract base class for detray detector payload output converters
class output_converter_interface {
 public:
  /// Default destructor
  virtual ~output_converter_interface() = default;

  /// Writes the detector payload to an output data sink.
  /// @param payload the detector payload to be added to the output data
  /// @param output the target output data
  virtual void from_payload(const detector_payload& payload,
                            const std::any& output) = 0;
};

/// @brief Abstract base class for detray detector payload file writers
class output_file_converter_interface : public output_converter_interface {
 public:
  /// Constructor setting the output format @param f and file extension
  /// @param ext
  output_file_converter_interface(const io::format f, const std::string& ext)
      : m_format{f}, m_file_extension{ext} {};

  /// @returns the output file format
  io::format format() const { return m_format; }

  /// @returns the file extension
  std::string_view file_extension() const { return m_file_extension; }

  /// @returns the mode in which the output files are opened
  std::ios_base::openmode open_mode() const { return m_mode; }

  /// Set the mode @param mode in which the output files are opened
  void open_mode(const std::ios_base::openmode mode) { m_mode = mode; }

  /// Writes the detector components in the payload to files.
  /// @note The existence of the file path has to be guaranteed by the caller
  ///
  /// @param payload the detector payload to be written to file
  /// @param path the path to the output file
  /// @param mode the open mode of the target file
  ///
  /// @returns the names of the files that were written
  virtual std::vector<std::string> write_to_file(
      const detector_payload& payload, const std::filesystem::path& path,
      const std::ios_base::openmode mode) = 0;

  /// Generic converter interface to the file writing
  /// @param payload the detector payload to be written to file
  /// @param output_data the path to the output file
  void from_payload(const detector_payload& payload,
                    const std::any& output_data) override {
    write_to_file(payload, get_file_path(output_data), m_mode);
  };

 private:
  /// @returns the output path that was passed as type-erased @param output_data
  std::filesystem::path get_file_path(const std::any& output_data) {
    std::filesystem::path file_path{"./"};
    try {
      file_path = std::any_cast<std::filesystem::path>(output_data);
    } catch (const std::bad_any_cast& e) {
      DETRAY_FATAL_HOST(
          "Unknown output data type in "
          << m_file_extension
          << " file converter! Should be 'std::filesystem::path'");
      throw;
    }

    return file_path;
  }

  /// Output file format of the respective writer
  io::format m_format;
  /// Extension that matches the file format of the respective writer
  std::string m_file_extension;
  /// How to open the output files
  std::ios_base::openmode m_mode{std::ios_base::out | std::ios_base::binary};
};

}  // namespace detray::io

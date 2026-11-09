// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/io/frontend/definitions.hpp"

// System include(s)
#include <ostream>
#include <string>

namespace detray::io {

/// @brief config struct for detector writing.
struct detector_writer_config {
  /// The path to the output files
  std::string m_path{"./"};
  /// The output file format
  detray::io::format m_format = detray::io::format::json;
  /// Replace files in case they already exists
  bool m_replace = false;
  /// Compactify json output, if not json format this flag does nothing
  bool m_compact_io = false;
  /// Whether to write the material to file
  bool m_write_material = true;
  /// Whether to write the accelerator grids to file
  bool m_write_grids = true;

  /// Getters
  /// @{
  /// Return the output directory used for detector files.
  /// @return output directory path
  const std::string& path() const { return m_path; }
  /// Return the selected output format.
  /// @return configured output format
  detray::io::format format() const { return m_format; }
  /// Return whether existing output files are replaced.
  /// @return whether replacement is enabled
  bool replace_files() const { return m_replace; }
  /// Return whether JSON output is compacted.
  /// @return whether JSON compaction is enabled
  bool compactify_json() const { return m_compact_io; }
  /// Return whether material data is written.
  /// @return whether material output is enabled
  bool write_material() const { return m_write_material; }
  /// Return whether accelerator grids are written.
  /// @return whether grid output is enabled
  bool write_grids() const { return m_write_grids; }
  /// @}

  /// Setters
  /// @{
  /// Set the output directory.
  /// @param p destination directory
  /// @return this configuration for chaining
  detector_writer_config& path(std::string p) {
    m_path = std::move(p);
    return *this;
  }
  /// Select the serialization format.
  /// @param f requested output format
  /// @return this configuration for chaining
  detector_writer_config& format(detray::io::format f) {
    m_format = f;
    return *this;
  }
  /// Configure whether to replace existing files.
  /// @param flag whether the option is enabled
  /// @return this configuration for chaining
  detector_writer_config& replace_files(bool flag) {
    m_replace = flag;
    return *this;
  }
  /// Configure whether to compactify JSON output.
  /// @param flag whether the option is enabled
  /// @return this configuration for chaining
  detector_writer_config& compactify_json(bool flag) {
    m_compact_io = flag;
    return *this;
  }
  /// Configure whether to include material data.
  /// @param flag whether the option is enabled
  /// @return this configuration for chaining
  detector_writer_config& write_material(bool flag) {
    m_write_material = flag;
    return *this;
  }
  /// Configure whether to include accelerator-grid data.
  /// @param flag whether the option is enabled
  /// @return this configuration for chaining
  detector_writer_config& write_grids(bool flag) {
    m_write_grids = flag;
    return *this;
  }
  /// @}

  /// Print the detector writer configuration
  friend std::ostream& operator<<(std::ostream& out,
                                  const detector_writer_config& cfg) {
    out << "\nDetector writer\n"
        << "----------------------------\n"
        << "  Path                  : " << cfg.path() << "\n"
        << "  Write grids           : " << std::boolalpha << cfg.write_grids()
        << "\n"
        << "  Write material        : " << cfg.write_material() << "\n";

    if (cfg.format() == detray::io::format::json) {
      out << "  Compactify json       : " << cfg.compactify_json() << "\n";
    }
    // Reset state
    out << std::noboolalpha;

    return out;
  }
};

}  // namespace detray::io

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

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
  const std::string& path() const { return m_path; }
  detray::io::format format() const { return m_format; }
  bool replace_files() const { return m_replace; }
  bool compactify_json() const { return m_compact_io; }
  bool write_material() const { return m_write_material; }
  bool write_grids() const { return m_write_grids; }
  /// @}

  /// Setters
  /// @{
  detector_writer_config& path(std::string p) {
    m_path = std::move(p);
    return *this;
  }
  detector_writer_config& format(detray::io::format f) {
    m_format = f;
    return *this;
  }
  detector_writer_config& replace_files(bool flag) {
    m_replace = flag;
    return *this;
  }
  detector_writer_config& compactify_json(bool flag) {
    m_compact_io = flag;
    return *this;
  }
  detector_writer_config& write_material(bool flag) {
    m_write_material = flag;
    return *this;
  }
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

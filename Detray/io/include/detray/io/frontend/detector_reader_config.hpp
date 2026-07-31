// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/volume_builder_interface.hpp"

// System include(s)
#include <initializer_list>
#include <ostream>
#include <string>
#include <vector>

namespace detray::io {

/// @brief config struct for detector reading.
struct detector_reader_config {
  /// Input files
  std::vector<std::string> m_files;
  /// Volume builder options
  detray::volume_builder_options m_builder_opts{};
  /// Run detector consistency check after reading
  bool m_do_check{true};
  /// Verbosity of the detector consistency check
  bool m_verbose{false};

  /// Getters
  /// @{
  const std::vector<std::string>& files() const { return m_files; }
  const detray::volume_builder_options& builder_options() const {
    return m_builder_opts;
  }
  bool deduplicate() const { return m_builder_opts.deduplicate(); }
  bool do_check() const { return m_do_check; }
  bool verbose_check() const { return m_verbose; }
  /// @}

  /// Setters
  /// @{
  detector_reader_config& add_file(const std::string& file_name) {
    m_files.push_back(file_name);
    return *this;
  }
  detector_reader_config& add_files(
      const std::vector<std::string>& file_names) {
    m_files = file_names;
    return *this;
  }
  detector_reader_config& add_files(std::vector<std::string>&& file_names) {
    m_files = std::move(file_names);
    return *this;
  }
  detector_reader_config& add_files(
      std::initializer_list<std::string> file_names) {
    m_files = std::vector<std::string>(file_names);
    return *this;
  }
  detector_reader_config& deduplicate(bool toggle) {
    m_builder_opts.deduplicate(toggle);
    return *this;
  }
  detector_reader_config& do_check(const bool check) {
    m_do_check = check;
    return *this;
  }
  detector_reader_config& verbose_check(const bool verbose) {
    if (verbose && !m_do_check) {
      m_do_check = true;
    }
    m_verbose = verbose;
    return *this;
  }
  /// @}

  /// Print the detector reader configuration
  friend std::ostream& operator<<(std::ostream& out,
                                  const detector_reader_config& cfg) {
    out << "\nDetector reader\n"
        << "----------------------------\n"
        << "  Detector files        : \n";
    for (const auto& file_name : cfg.files()) {
      out << "    -> " << file_name << "\n";
    }
    out << "  Deduplicate data      : " << std::boolalpha << cfg.deduplicate()
        << std::noboolalpha << "\n";

    return out;
  }
};

}  // namespace detray::io

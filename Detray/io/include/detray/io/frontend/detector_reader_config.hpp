// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

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
  /// Run detector consistency check after reading
  bool m_do_check{true};
  /// Verbosity of the detector consistency check
  bool m_verbose{false};

  /// Getters
  /// @{
  const std::vector<std::string>& files() const { return m_files; }
  bool do_check() const { return m_do_check; }
  bool verbose_check() const { return m_verbose; }
  /// @}

  /// Setters
  /// @{
  detector_reader_config& add_file(std::string_view file_name) {
    m_files.push_back(std::string{file_name});
    return *this;
  }
  detector_reader_config& add_files(
      const std::vector<std::string>& file_names) {
    m_files.reserve(m_files.size() + file_names.size());
    m_files.insert(std::end(m_files), std::begin(file_names),
                   std::end(file_names));
    return *this;
  }
  detector_reader_config& add_files(
      const std::vector<std::string_view>& file_names) {
    m_files.reserve(m_files.size() + file_names.size());
    m_files.insert(std::end(m_files), std::begin(file_names),
                   std::end(file_names));
    return *this;
  }
  detector_reader_config& add_files(std::vector<std::string>&& file_names) {
    m_files.reserve(m_files.size() + file_names.size());
    m_files.insert(std::end(m_files),
                   std::make_move_iterator(std::begin(file_names)),
                   std::make_move_iterator(std::end(file_names)));
    return *this;
  }
  template <typename... Args>
    requires(std::same_as<std::string, std::remove_cvref_t<Args>> || ...)
  detector_reader_config& add_files(Args&&... file_names) {
    std::array<std::string, sizeof...(Args)> tmp_files{file_names...};

    m_files.reserve(m_files.size() + tmp_files.size());
    m_files.insert(std::end(m_files), std::begin(tmp_files),
                   std::end(tmp_files));
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
        << "  Detector files:       : \n";
    for (const auto& file_name : cfg.files()) {
      out << "    -> " << file_name << "\n";
    }

    return out;
  }
};

}  // namespace detray::io

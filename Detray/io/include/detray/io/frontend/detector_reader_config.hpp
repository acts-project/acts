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
  /// Return the list of configured detector input files.
  /// @return configured file paths
  const std::vector<std::string>& files() const { return m_files; }
  /// Return the volume-builder options used during reconstruction.
  /// @return configured builder options
  const detray::volume_builder_options& builder_options() const {
    return m_builder_opts;
  }
  /// Return whether duplicate detector components are removed.
  /// @return whether the option is enabled
  bool deduplicate() const { return m_builder_opts.deduplicate(); }
  /// Return whether detector consistency checks are enabled.
  /// @return whether the option is enabled
  bool do_check() const { return m_do_check; }
  /// Return whether consistency checks print detailed output.
  /// @return whether the option is enabled
  bool verbose_check() const { return m_verbose; }
  /// @}

  /// Setters
  /// @{
  /// Add a detector component file to the input list.
  /// @param file_name path to a detector data file
  /// @return this configuration for chaining
  detector_reader_config& add_file(const std::string& file_name) {
    m_files.push_back(file_name);
    return *this;
  }
  /// Enable or disable deduplication in the volume builder.
  /// @param toggle whether duplicate components should be removed
  /// @return this configuration for chaining
  detector_reader_config& deduplicate(bool toggle) {
    m_builder_opts.deduplicate(toggle);
    return *this;
  }
  /// Enable or disable the post-read detector consistency check.
  /// @param check whether to run the check
  /// @return this configuration for chaining
  detector_reader_config& do_check(const bool check) {
    m_do_check = check;
    return *this;
  }
  /// Enable detailed output from the consistency check.
  /// Enabling verbose output also enables the check.
  /// @param verbose whether to print detailed check output
  /// @return this configuration for chaining
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

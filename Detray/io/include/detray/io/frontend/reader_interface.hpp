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

// System include(s)
#include <filesystem>
#include <ios>
#include <string>
#include <string_view>

namespace detray::io {

/// @brief Abstract base class for detray detector component readers
template <class detector_t>
class reader_interface {
 public:
  /// All readers must define a file extension
  reader_interface() = delete;

  /// Only accept files with a fixed @param extension
  explicit reader_interface(const std::string& ext) : m_file_extension{ext} {}

  /// Default destructor
  virtual ~reader_interface() = default;

  /// @returns the file extension
  const std::string& file_extension() const { return m_file_extension; }

  /// Reads the respective detector component from file. Since the detector
  /// does not keep the volume names, the name map is also passed and
  /// filled.
  virtual void read(
      detector_builder<typename detector_t::metadata, volume_builder>&,
      const std::string&) = 0;

 private:
  /// Extension that matches the file format of the respective reader
  std::string m_file_extension;
};

}  // namespace detray::io

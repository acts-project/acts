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

/// @brief Abstract base class for detray detector component writers
template <class detector_t>
class writer_interface {
 public:
  /// All writers must define a file extension
  writer_interface() = delete;

  /// File gets created with a fixed @param extension
  explicit writer_interface(const std::string& ext) : m_file_extension{ext} {}

  /// Default destructor
  virtual ~writer_interface() = default;

  /// @returns the file extension
  const std::string& file_extension() const { return m_file_extension; }

  /// Writes the respective detector component to file. Since the detector
  /// does not provide the volume names, the name map is also passed.
  /// @note The existence of the file path has to be guaranteed by the caller
  virtual std::string write(const detector_t&,
                            const typename detector_t::name_map&,
                            const std::ios_base::openmode,
                            const std::filesystem::path&) = 0;

 private:
  /// Extension that matches the file format of the respective writer
  std::string m_file_extension;
};

}  // namespace detray::io

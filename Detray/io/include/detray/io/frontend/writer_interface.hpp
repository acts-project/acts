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

  /// Write files with the given extension.
  /// @param ext file extension handled by this writer
  explicit writer_interface(const std::string& ext) : m_file_extension{ext} {}

  /// Default destructor
  virtual ~writer_interface() = default;

  /// @returns the file extension
  const std::string& file_extension() const { return m_file_extension; }

  /// Write detector components to a path using the supplied name map.
  /// @param detector detector components to serialize
  /// @param names map used to name detector components
  /// @param mode stream open mode
  /// @param path output directory guaranteed to exist by the caller
  /// @return path or filename written by the implementation
  virtual std::string write(const detector_t& detector,
                            const typename detector_t::name_map& names,
                            const std::ios_base::openmode mode,
                            const std::filesystem::path& path) = 0;

 private:
  /// Extension that matches the file format of the respective writer
  std::string m_file_extension;
};

}  // namespace detray::io

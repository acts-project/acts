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
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/io/frontend/reader_interface.hpp"

// System include(s)
#include <algorithm>
#include <cassert>
#include <map>
#include <memory>
#include <string>
#include <type_traits>

namespace detray::io::detail {

/// @brief A reader for multiple detector components.
///
/// The class aggregates a number of different readers and calls them once the
/// detector data should be read in from file.
template <class detector_t>
class detector_components_reader final {
  using reader_ptr_t = std::unique_ptr<reader_interface<detector_t>>;

 public:
  /// Default constructor
  detector_components_reader() = default;

  /// Create a new reader of type @tparam reader_t
  template <class reader_t>
    requires std::is_base_of_v<reader_interface<detector_t>, reader_t>
  void add(const std::string& file_name) {
    add(std::make_unique<reader_t>(), file_name);
  }

  /// Attach an existing reader via @param r_ptr to the readers
  void add(reader_ptr_t&& r_ptr, const std::string& file_name) {
    m_readers[file_name] = std::move(r_ptr);
  }

  /// Create a new geometry reader of type @tparam reader_t
  template <class reader_t>
  void add_geometry_reader(std::string_view file_name) {
    m_geo_file_name = file_name;
    m_geo_reader = std::make_unique<reader_t>();
  }

  /// @returns the number of readers that are registered
  std::size_t size() const {
    return m_readers.size() + (m_geo_reader == nullptr ? 0u : 1u);
  }

  /// @returns access to the readers map - const
  const auto& readers_map() const { return m_readers; }

  /// Set the name of the detector to be read
  void set_detector_name(std::string name) { m_det_name = std::move(name); }

  /// Reads the full detector into @param det by calling the readers, while
  /// using the name map @param volume_names for to write the volume names.
  void read(detector_builder<typename detector_t::metadata, volume_builder>&
                det_builder) {
    // We have to at least read a geometry
    assert(size() != 0u &&
           "No readers registered! Need at least a geometry reader");

    // Set the detector name in the name map
    det_builder.set_name(m_det_name);

    if (m_geo_reader == nullptr) {
      throw std::runtime_error("No geometry reader registered!");
    }
    m_geo_reader->read(det_builder, m_geo_file_name);

    // Call the read method on all readers
    for (const auto& [file_name, reader] : m_readers) {
      reader->read(det_builder, file_name);
    }
  }

 private:
  /// Name of the detector
  std::string m_det_name;
  /// Geometry file name
  std::string m_geo_file_name;
  /// Geometry reader
  reader_ptr_t m_geo_reader;
  /// The readers registered for the detector: e.g. material, grids...
  std::unordered_map<std::string, reader_ptr_t> m_readers;
};

}  // namespace detray::io::detail

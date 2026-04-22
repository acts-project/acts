// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/io/frontend/writer_interface.hpp"

// System include(s)
#include <algorithm>
#include <cassert>
#include <ios>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace detray::io::detail {

/// @brief A writer for multiple detector components.
///
/// The class aggregates a number of different writers and calls them once the
/// detector data should be written to file. The respective writers are
/// responsible for file handling etc.
///
/// @note the writers are unordered and must not depend on a write order!
template <class detector_t>
class detector_components_writer final {
  using writer_ptr_t = std::unique_ptr<writer_interface<detector_t>>;

 public:
  /// Default constructor
  detector_components_writer() = default;

  /// Create a new writer of type @tparam writer_t
  template <class writer_t>
    requires std::is_base_of_v<writer_interface<detector_t>, writer_t>
  void add() {
    add(std::make_unique<writer_t>());
  }

  /// Attach an existing writer via @param w_ptr to the writers
  void add(writer_ptr_t&& w_ptr) { m_writers.push_back(std::move(w_ptr)); }

  /// Writes the full detector data of @param det to file by calling the
  /// writers, while using the name map @param names for the detector
  void write(const detector_t& det, const typename detector_t::name_map& names,
             const std::ios_base::openmode mode,
             const std::filesystem::path& file_path) {
    // We have to at least write a geometry
    assert(!m_writers.empty() &&
           "No writers registered! Need at least a geometry writer");

    // Call the write method on all optional writers
    std::ranges::for_each(
        m_writers, [&det, &names, mode, &file_path](writer_ptr_t& writer) {
          writer->write(det, names, mode, file_path);
        });
  }

 private:
  /// The writers registered for the detector: geometry (mandatory!) plus
  /// e.g. material, grids...)
  std::vector<writer_ptr_t> m_writers;
};

}  // namespace detray::io::detail

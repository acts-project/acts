// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/io/read_detector.hpp"

#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/detector.hpp"
#include "traccc/io/utils.hpp"

// Detray include(s)
#include <detray/core/concepts.hpp>

// System include(s).
#include <string>
#include <string_view>
#include <tuple>

namespace {

/// Common implementation for constructing a detector from a set of input files
template <detray::concepts::detector detector_t>
void read_detector(traccc::host_detector& detector, vecmem::memory_resource& mr,
                   const std::string_view& geometry_file,
                   const std::string_view& material_file,
                   const std::string_view& grid_file,
                   const bool do_consistency_check = false) {
  // Set up the detector reader configuration.
  detray::io::detector_reader_config cfg;
  cfg.do_check(do_consistency_check);

  cfg.add_file(traccc::io::get_absolute_path(geometry_file));
  if (material_file.empty() == false) {
    cfg.add_file(traccc::io::get_absolute_path(material_file));
  }
  if (grid_file.empty() == false) {
    cfg.add_file(traccc::io::get_absolute_path(grid_file));
  }

  // Read the detector.
  auto det = detray::io::read_detector<detector_t>(mr, cfg);
  detector.set<detector_t>(std::move(det.first));
}

/// Check if a detector name matches "<name>" or "<name>_detector"
///
/// @param detector_name the detector name from the geometry file header
/// @param name the name of a generated detector type
///
/// @returns @c true if the names match
[[maybe_unused]] bool matches_detector_name(
    const std::string_view detector_name, const std::string_view name) {
  return (detector_name == name) ||
         (detector_name == std::string(name) + "_detector");
}

/// Construct the generated detector type that matches a detector name
///
/// @tparam detector_ts the generated detector types to check
///
/// @param detector_name the detector name from the geometry file header
/// @param args the arguments for @c read_detector
///
/// @returns @c true if a matching detector type was found and constructed
/// (always @c false if no detector type is generated)
template <typename... detector_ts, typename... arg_ts>
bool read_generated_detector(
    std::tuple<detector_ts...>*,
    [[maybe_unused]] const std::string_view detector_name,
    [[maybe_unused]] arg_ts&&... args) {
  // Read the first matching detector type
  return ((matches_detector_name(
               detector_name, traccc::generated_detector_name<detector_ts>) &&
           (::read_detector<detector_ts>(args...), true)) ||
          ...);
}

}  // namespace

namespace traccc::io {

void read_detector(host_detector& detector, vecmem::memory_resource& mr,
                   const std::string_view& geometry_file,
                   const std::string_view& material_file,
                   const std::string_view& grid_file) {
  // Peek at the header to determine the kind of detector that is needed
  const auto header = detray::io::detail::deserialize_json_header(
      traccc::io::get_absolute_path(geometry_file));

  // TODO: Update this
  if (header.detector == "Cylindrical detector from DD4hep blueprint") {
    ::read_detector<traccc::odd_detector>(detector, mr, geometry_file,
                                          material_file, grid_file);
  } else if (!::read_generated_detector(
                 static_cast<generated_detector_type_list*>(nullptr),
                 header.detector, detector, mr, geometry_file, material_file,
                 grid_file)) {
    // TODO: Warning here
    ::read_detector<traccc::default_detector>(detector, mr, geometry_file,
                                              material_file, grid_file);
  }
}

}  // namespace traccc::io

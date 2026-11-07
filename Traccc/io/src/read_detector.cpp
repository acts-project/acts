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
#include <detray/utils/type_list.hpp>

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
                   const std::string_view& grid_file) {
  // Set up the detector reader configuration.
  detray::io::detector_reader_config cfg;
  cfg.do_check(false);

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

/// Construct the generated detector type that matches a metadata type
///
/// @tparam detector_ts the generated detector types to check
///
/// @param metadata_name the metadata type from the geometry file header
/// @param args the arguments for @c read_detector
///
/// @returns @c true if a matching detector type was found and constructed
/// (always @c false if no detector type is generated)
template <typename... metadata_ts, typename... arg_ts>
bool read_generated_detector(
    std::tuple<metadata_ts...>*,
    [[maybe_unused]] const std::string_view metadata_name,
    [[maybe_unused]] arg_ts&&... args) {
  // Read the first matching detector type
  return ((metadata_name == detray::types::get_name<metadata_ts>() &&
           (::read_detector<detray::host::detector<metadata_ts>>(
                std::forward<arg_ts>(args)...),
            true)) ||
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

  // TODO: Update this and remove detector name
  const std::string_view det_name{header.detector};
  const std::string_view metadata{header.metadata};

  // TODO: Remove special cases, once the json files are fixed
  if (metadata == "odd_metadata" ||
      det_name == "Cylindrical detector from DD4hep blueprint") {
    ::read_detector<traccc::odd_detector>(detector, mr, geometry_file,
                                          material_file, grid_file);
#ifdef DETRAY_HAS_ITK_METADATA
  } else if (metadata == "itk_metadata" || det_name == "detray_detector") {
    ::read_detector<traccc::itk_detector>(detector, mr, geometry_file,
                                          material_file, grid_file);
#endif
  } else if (!::read_generated_detector(
                 static_cast<generated_metadata_list*>(nullptr),
                 header.metadata, detector, mr, geometry_file, material_file,
                 grid_file)) {
    // TODO: Warning here
    ::read_detector<traccc::default_detector>(detector, mr, geometry_file,
                                              material_file, grid_file);
  }
}

}  // namespace traccc::io

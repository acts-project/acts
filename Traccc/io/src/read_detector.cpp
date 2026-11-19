// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/io/read_detector.hpp"

#include <detray/detectors/io.hpp>

#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/detector.hpp"
#include "traccc/io/utils.hpp"

namespace traccc::io {

void read_detector(host_detector& detector, vecmem::memory_resource& mr,
                   const std::string_view geometry_file,
                   const std::string_view material_file,
                   const std::string_view grid_file,
                   const bool do_consistency_check,
                   const std::string_view metadata_id) {
  detray::io::detector_reader_config cfg;
  cfg.do_check(do_consistency_check)
      .add_files(
          std::vector<std::string>{traccc::io::get_absolute_path(geometry_file),
                                   traccc::io::get_absolute_path(material_file),
                                   traccc::io::get_absolute_path(grid_file)});
  auto result = detray::io::read_registered_detector_json<default_algebra>(
      mr, cfg, metadata_id);
  detector = std::move(result.first);
}

}  // namespace traccc::io

// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Acts include(s)
#include "Acts/Geometry/GeometryIdentifier.hpp"

// Acts Examples include(s)
#include "ActsExamples/EventData/Cluster.hpp"

#include <detray/io/frontend/detector_reader.hpp>

namespace ActsExamples::Traccc::Common::Io {
template <typename detector_t>
inline auto readDetector(vecmem::memory_resource* mr,
                         const std::string& detectorFilePath,
                         const std::string& materialFilePath = "",
                         const std::string& gridFilePath = "") {
  // Set up the detector reader configuration.
  detray::io::detector_reader_config cfg;
  cfg.add_file(detectorFilePath);
  if (!materialFilePath.empty()) {
    cfg.add_file(materialFilePath);
  }
  if (!gridFilePath.empty()) {
    cfg.add_file(gridFilePath);
  }

  // Read the detector.
  auto [det, names] = detray::io::read_detector<detector_t>(*mr, cfg);
  return std::move(det);
}
}  // namespace ActsExamples::Traccc::Common::Io

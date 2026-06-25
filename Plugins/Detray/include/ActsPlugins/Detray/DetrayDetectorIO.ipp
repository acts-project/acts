// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Implementation of the detray detector operations. Include this header
// (instead of DetrayDetectorIO.hpp) when instantiating the operations for a
// metadata type that is not part of the closed set declared in
// DetrayMetadata.hpp.

#include "ActsPlugins/Detray/DetrayDetectorIO.hpp"

#include <utility>

#include <detray/io/frontend/detector_reader.hpp>
#include <detray/io/frontend/detector_reader_config.hpp>
#include <detray/io/frontend/detector_writer.hpp>
#include <detray/io/frontend/detector_writer_config.hpp>
#include <detray/utils/consistency_checker.hpp>

namespace ActsPlugins::detail {

template <typename metadata_t>
void checkDetrayConsistency(const detray::detector<metadata_t>& detector) {
  detray::detail::check_consistency(detector);
}

template <typename metadata_t>
void writeDetrayJson(const detray::detector<metadata_t>& detector,
                     const detray::name_map& names, const std::string& path) {
  auto cfg = detray::io::detector_writer_config{}
                 .format(detray::io::format::json)
                 .path(path)
                 .replace_files(true);
  detray::io::write_detector(detector, names, cfg);
}

template <typename metadata_t>
std::pair<detray::detector<metadata_t>, detray::name_map> readDetrayDetector(
    vecmem::memory_resource& mr, const std::vector<std::string>& files) {
  auto cfg = detray::io::detector_reader_config{}.do_check(false);
  for (const auto& file : files) {
    cfg.add_file(file);
  }
  return detray::io::read_detector<detray::detector<metadata_t>>(mr, cfg);
}

}  // namespace ActsPlugins::detail

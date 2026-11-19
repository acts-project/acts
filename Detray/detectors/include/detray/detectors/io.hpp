// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "detray/core/host_detector.hpp"
#include "detray/detectors/metadata.hpp"
#include "detray/io/frontend/detector_writer.hpp"
#include "detray/io/json/detector_reader.hpp"

namespace detray::io {

/// Construct a registered detector from a payload. An explicit ID overrides the
/// payload header, allowing callers to read legacy files with missing metadata.
/// Unknown/missing IDs never fall back to a different detector family.
/// The supplied memory resource must outlive the returned detector.
template <typename algebra_t>
std::pair<host_detector, name_map> read_registered_detector(
    vecmem::memory_resource& mr, const detector_reader_config& cfg,
    detector_payload& payload, std::string_view metadata_id = {}) {
  const auto id = metadata_id.empty()
                      ? std::string_view{payload.header.metadata}
                      : metadata_id;
  return registry::visit_metadata<algebra_t>(
      id,
      [&]<detray::concepts::metadata metadata_t>()
          -> std::pair<host_detector, name_map> {
        auto [concrete, names] =
            read_detector<host::detector<metadata_t>>(mr, cfg, payload);
        host_detector detector;
        detector.set<metadata_t>(std::move(concrete));
        return {std::move(detector), std::move(names)};
      });
}

/// Read JSON files using their registered metadata ID, or an explicit override.
template <typename algebra_t>
std::pair<host_detector, name_map> read_registered_detector_json(
    vecmem::memory_resource& mr, const detector_reader_config& cfg,
    std::string_view metadata_id = {}) {
  detector_payload payload;
  convert_to_payload<json_input_converter>(payload, cfg.files());
  auto build_config = cfg;
  build_config.m_files
      .clear();  // Files have already been converted to payloads.
  return read_registered_detector<algebra_t>(mr, build_config, payload,
                                             metadata_id);
}

/// Write a registered detector with its stable family ID in every component.
template <typename algebra_t>
void write_registered_detector(const host_detector& detector,
                               const name_map& names,
                               detector_writer_config& cfg) {
  host_detector_visitor<registry::detector_types<algebra_t>>(
      detector,
      [&]<detray::concepts::detector detector_t>(const detector_t& concrete) {
        write_detector(
            concrete, names, cfg,
            registry::metadata_id<detector_metadata_t<detector_t>>());
      });
}

}  // namespace detray::io

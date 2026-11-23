// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/io/frontend/detector_writer.hpp"
#include "detray/io/json/json_converter.hpp"

namespace detray::io {

/// @brief Convert the detector payload to json output files
///
/// @param[in, out] payload payload containing detector data (all components).
///                         Its header is generated from the payload data
/// @param[in] cfg the detector writer configuration
inline void convert_payload_to_json(
    detray::io::detector_payload& payload,
    const detray::io::detector_writer_config& cfg) {
  detray::io::convert_from_payload<json_output_converter>(payload, cfg);
}

/// @brief Write the detector completely to json files.
///
/// @tparam detector_t the type of detector to be written
///
/// @param det the detector to be written
/// @param names the volume names of the detector
/// @param cfg the detector writer configuration
template <class detector_t>
void write_detector_json(
    const detector_t& det, const typename detector_t::name_map& names,
    detray::io::detector_writer_config& cfg) noexcept(false) {
  detray::io::write_detector<detector_t, json_output_converter>(det, names,
                                                                cfg);
}

}  // namespace detray::io

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/io/json/json_converter.hpp"

namespace detray::io {

/// @brief Convert data from json input files to detector payload
///
/// @param[out] payload payload containing detector data (all components)
/// @param[in] file_names a number of input file names
void convert_json_to_payload(detray::io::detector_payload& payload,
                             std::initializer_list<std::string> file_names) {
  detray::io::convert_to_payload<json_input_converter>(payload, {file_names});
}

/// @brief Read and build the detector completely from json files.
///
/// @tparam detector_t the type of detector to be built
/// @tparam GCAP surface grid bin capacity. If CAP is 0, the grid reader builds
///              a grid type with dynamic bin capacity
/// @tparam GDIM dimension of the surface grids, usually 2D
/// @tparam MDIM dimension of the material grids, usually 2D
/// @tparam volume_builder_t the type of base volume builder to be used
///
/// @param resc the memory resource to be used for the detector container allocs
/// @param cfg the detector reader configuration
///
/// @returns a complete detector object + a map that contains the volume names
template <class detector_t, std::size_t GCAP = 0u, std::size_t GDIM = 2u,
          std::size_t MDIM = 2u,
          template <typename> class volume_builder_t = volume_builder>
auto read_detector_json(vecmem::memory_resource& resc,
                        const detector_reader_config& cfg) noexcept(false) {
  return read_detector<detector_t, GCAP, GDIM, MDIM, volume_builder_t,
                       json_input_converter>(resc, cfg);
}

}  // namespace detray::io

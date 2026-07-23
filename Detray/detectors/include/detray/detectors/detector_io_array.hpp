// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Extern-template facade for the array-algebra detray detectors.
//
// Including this header opts every shipped `detector<meta<array<float>>>` into
// the pre-built IO operations: the heavy `write_detector`, `read_detector`,
// `check_consistency` and `assemble_detector` template trees are compiled once
// in the `detray::detector_io_array` library and merely linked here, instead of
// being instantiated in every consumer.
//
// This is the array counterpart of the `detray::io` vs `detray::io_array`
// split: the frontend/core headers cannot carry these declarations themselves
// because they live below the `detectors` layer and cannot name the concrete
// metadata types. Consumers that need a non-shipped metadata or a different
// algebra simply include the frontend headers directly and instantiate on use.

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/detectors/default_metadata.hpp"
#include "detray/detectors/detail/detector_io_instantiations.hpp"
#include "detray/detectors/odd_metadata.hpp"
#include "detray/detectors/telescope_metadata.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/detectors/wire_chamber_metadata.hpp"
#include "detray/io/frontend/definitions.hpp"
#include "detray/io/frontend/detector_assembler.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/io/frontend/detector_writer.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/utils/consistency_checker.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <string_view>
#include <utility>

// clang-format off
#define DETRAY_DETECTOR_IO_EXTERN(META)                                        \
  extern template void detray::io::write_detector<                             \
      detray::detector<detray::META<detray::array<float>>>>(                   \
      detray::detector<detray::META<detray::array<float>>>&,                   \
      const detray::detector<detray::META<detray::array<float>>>::name_map&,   \
      detray::io::detector_writer_config&);                                    \
  extern template std::pair<                                                   \
      detray::detector<detray::META<detray::array<float>>>,                    \
      detray::detector<detray::META<detray::array<float>>>::name_map>          \
  detray::io::read_detector<                                                   \
      detray::detector<detray::META<detray::array<float>>>>(                   \
      vecmem::memory_resource&, const detray::io::detector_reader_config&);    \
  extern template bool detray::detail::check_consistency<                      \
      detray::detector<detray::META<detray::array<float>>>>(                   \
      const detray::detector<detray::META<detray::array<float>>>&, bool,       \
      const detray::detector<detray::META<detray::array<float>>>::name_map&);  \
  extern template detray::detector<detray::META<detray::array<float>>>         \
  detray::io::assemble_detector<                                               \
      detray::detector<detray::META<detray::array<float>>>>(                   \
      vecmem::memory_resource&, const detray::io::detector_payload&,           \
      const detray::io::detector_homogeneous_material_payload*,                \
      detray::io::detector_grids_payload<                                      \
          detray::io::surface_material_payload, detray::io::material_id>*,     \
      const detray::io::detector_grids_payload<std::size_t,                    \
                                               detray::io::accel_id>*,         \
      std::string_view,                                                        \
      detray::detector<detray::META<detray::array<float>>>::name_map&);
// clang-format on

DETRAY_IO_METADATA_FOR_EACH(DETRAY_DETECTOR_IO_EXTERN)

#undef DETRAY_DETECTOR_IO_EXTERN

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"

// Detray include(s).

// Turn off excessive detector building debug logs (only works with gcc!)
// clang-format off
#include <detray/utils/quiet_log_start.hpp>
#include <detray/core/detector.hpp>
#include <detray/utils/quiet_log_end.hpp>
// clang-format on

#include <detray/detectors/default_metadata.hpp>
// Metadata that was generated in this build (see DETRAY_GENERATE_METADATA)
#include <detray/detectors/generated_metadata.hpp>
#include <detray/detectors/odd_metadata.hpp>
#include <detray/detectors/telescope_metadata.hpp>
#include <detray/detectors/toy_metadata.hpp>
#include <detray/detectors/wire_chamber_metadata.hpp>
#include <detray/utils/tuple_helpers.hpp>

// System include(s).
#include <string_view>
#include <tuple>

namespace traccc {

/// Default detector (can contain (almost) any detector data)
using default_detector =
    detray::host::detector<detray::default_metadata<traccc::default_algebra>>;

using default_detector_traits = detray::detector_traits<default_detector>;

/// Open Data Detector (ODD) detector
using odd_detector =
    detray::host::detector<detray::odd_metadata<traccc::default_algebra>>;

/// Detray telescope detector (test detector)
using telescope_detector = detray::host::detector<
    detray::telescope_metadata<traccc::default_algebra, detray::rectangle2D>>;

/// Detray toy detector (test detector)
using toy_detector =
    detray::host::detector<detray::toy_metadata<traccc::default_algebra>>;

/// Detray wire chamber detector (test detector)
using wire_chamber = detray::host::detector<
    detray::wire_chamber_metadata<traccc::default_algebra>>;

/// Name of a generated detector type
///
/// For the metadata @c detray::<name>_metadata the name is @c "<name>". It is
/// used to match the detector name in a geometry file header.
///
/// @tparam detector_t the detector type
template <typename detector_t>
inline constexpr std::string_view generated_detector_name{};

/// Detector types for all metadata generated in this build
///
/// For each metadata @c detray::<name>_metadata (see
/// @c DETRAY_GENERATE_METADATA), define the detector type @c <name>_detector
/// and its name @c generated_detector_name.
///
/// @note The generated metadata must be templated on the algebra type.
/// @{
// clang-format off
#define TRACCC_DEFINE_GENERATED_DETECTOR(NAME)                          \
  using NAME##_detector = detray::host::detector<                       \
      detray::NAME##_metadata<traccc::default_algebra>>;                \
  template <>                                                           \
  inline constexpr std::string_view                                     \
      generated_detector_name<NAME##_detector> = #NAME;
DETRAY_FOR_EACH_GENERATED_METADATA(TRACCC_DEFINE_GENERATED_DETECTOR)
#undef TRACCC_DEFINE_GENERATED_DETECTOR

#define TRACCC_GENERATED_DETECTOR_TUPLE(NAME) , std::tuple<NAME##_detector>
/// List of all generated detector types (can be empty)
using generated_detector_type_list = detray::detail::tuple_cat_t<
    std::tuple<> DETRAY_FOR_EACH_GENERATED_METADATA(
        TRACCC_GENERATED_DETECTOR_TUPLE)>;
#undef TRACCC_GENERATED_DETECTOR_TUPLE

#define TRACCC_GENERATED_METADATA_TUPLE(NAME) , std::tuple<detray::NAME##_metadata<traccc::default_algebra>>
/// List of all generated detector types (can be empty)
using generated_metadata_list = detray::detail::tuple_cat_t<
    std::tuple<> DETRAY_FOR_EACH_GENERATED_METADATA(
        TRACCC_GENERATED_METADATA_TUPLE)>;
#undef TRACCC_GENERATED_METADATA_TUPLE
// clang-format on
/// @}

}  // namespace traccc

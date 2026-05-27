/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"

// Detray include(s).

// Turn off exessive detector building debug logs (only works with gcc!)
// clang-format off
#include <detray/utils/quiet_log_start.hpp>
#include <detray/core/detector.hpp>
#include <detray/utils/quiet_log_end.hpp>
// clang-format on

#include <detray/detectors/default_metadata.hpp>
#include <detray/detectors/itk_metadata.hpp>
#include <detray/detectors/odd_metadata.hpp>
#include <detray/detectors/telescope_metadata.hpp>
#include <detray/detectors/toy_metadata.hpp>
#include <detray/detectors/wire_chamber_metadata.hpp>

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

// System include(s).
#include <type_traits>

namespace traccc {
namespace details {

/// Type used instead of @c detray::device_container_types
///
/// This is meant as a template for how @c detray::device_container_types should
/// change. So that it would correctly reflect that all data access on a device
/// is constant. Not modifying the detector's payload.
///
/// Also note that with all the types present in @c detray::container_types
/// it's really only @c vector_type that is actually used by @c detray::detector
/// at this point.
///
struct device_detector_container_types {

    /// Vector type to use in device code
    template <typename T>
    using vector_type = vecmem::device_vector<std::add_const_t<T>>;

};  // struct device_detector_container_types

}  // namespace details

/// Base struct for the different detector types supported by the project.
template <typename metadata_t>
struct detector_traits {

    /// Metadata type of the detector.
    using metadata_type = metadata_t;

    /// Host type of the detector.
    using host = detray::detector<metadata_type, detray::host_container_types>;
    /// Device type of the detector.
    using device = detray::detector<metadata_type,
                                    details::device_detector_container_types>;

    /// Non-const view of the detector.
    using view = typename host::const_view_type;

    /// Buffer for a detector's data.
    using buffer = typename host::buffer_type;

};  // struct default_detector

template <typename T>
concept is_detector_traits = requires {
    typename T::metadata_type;
    typename T::host;
    typename T::device;
    typename T::view;
    typename T::buffer;
};

/// Default detector (can contain any detector data)
using default_detector =
    detector_traits<detray::default_metadata<traccc::default_algebra>>;

/// ATLAS Inner Tracker (ITk) detector
using itk_detector =
    detector_traits<detray::itk_metadata<traccc::default_algebra>>;

/// Open Data Detector (ODD) detector
using odd_detector =
    detector_traits<detray::odd_metadata<traccc::default_algebra>>;

/// Detray telescope detector (test detector)
using telescope_detector = detector_traits<
    detray::telescope_metadata<traccc::default_algebra, detray::rectangle2D>>;

/// Detray toy detector (test detector)
using toy_detector =
    detector_traits<detray::toy_metadata<traccc::default_algebra>>;

/// Detray wire chamber detector (test detector)
using wire_chamber =
    detector_traits<detray::wire_chamber_metadata<traccc::default_algebra>>;

}  // namespace traccc

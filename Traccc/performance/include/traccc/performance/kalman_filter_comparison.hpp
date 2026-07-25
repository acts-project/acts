/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/particle.hpp"
#include "traccc/utils/seed_generator.hpp"
#include "traccc/utils/transcribe_to_trace.hpp"

// Detray include(s)
#include <detray/propagator/propagation_config.hpp>
#include <detray/test/validation/propagation_validation.hpp>

// System include(s)
#include <memory>
#include <string>

namespace traccc {

/// Run a detray propagation with Kalman filter and compare the surfaces
/// and measurements that were encountered against a reference of
/// truth tracks
///
/// @param det the detector
/// @param names the detector and volume names
/// @param bfield the magnetic field representation
/// @param cfg the full test configuration
/// @param ilogger logging service
/// @param tracks the initial particle vertices as free track parameters
/// @param truth_traces_fw the reference data for each encountered module
/// @param device_measurements the truth measurements for the Kalman Filter
/// @param track_container the truth track container for the Kalman Filter
///
/// @returns whether the validation was successful
bool kalman_filter_comparison(
    const traccc::default_detector::host& det,
    const traccc::default_detector::host::name_map& names,
    const traccc::magnetic_field& bfield,
    const detray::propagation_validation_config& cfg,
    const traccc::seed_generator<traccc::default_detector::host>::config&
        smearing_cfg,
    std::unique_ptr<const traccc::Logger> ilogger,
    const std::vector<traccc::free_track_parameters<traccc::default_algebra>>&
        tracks,
    std::vector<vecmem::vector<traccc::propagation_validator::candidate_type<
        traccc::default_detector::host>>>& truth_traces_fw,
    edm::measurement_collection::const_device device_measurements,
    traccc::edm::track_container<traccc::default_algebra>::host&
        track_container);

}  // namespace traccc

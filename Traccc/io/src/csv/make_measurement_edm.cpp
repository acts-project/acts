/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/csv/make_measurement_edm.hpp"

// Project include(s).
#include "traccc/definitions/track_parametrization.hpp"

// Detray include(s).
#include <detray/geometry/identifier.hpp>

namespace traccc::io::csv {

void make_measurement_edm(
    const traccc::io::csv::measurement& csv_meas,
    edm::measurement_collection::host::proxy_type& meas,
    const std::map<geometry_id, geometry_id>* acts_to_detray_id,
    const traccc::detector_design_description::host* det_desc,
    const std::map<geometry_id, std::size_t>*
        geometry_id_to_detector_description_index) {

    // Construct the measurement object.
    std::array<detray::dindex_type<default_algebra>, 2u> indices{0u, 0u};
    meas.dimensions() = 0u;

    // Local key is a 8 bit char and first and last bit are dummy value. 2 -
    // 7th bits are for 6 bound track parameters.
    // Ex1) 0000010 or 2 -> meas dim = 1 and [loc0] active -> strip or wire
    // Ex2) 0000110 or 6 -> meas dim = 2 and [loc0, loc1] active -> pixel
    // Ex3) 0000100 or 4 -> meas dim = 1 and [loc1] active -> annulus
    for (unsigned int ipar = 0; ipar < 2u; ++ipar) {
        if (((csv_meas.local_key) & (1 << (ipar + 1))) != 0) {

            switch (ipar) {
                case e_bound_loc0: {
                    meas.local_position()[0] = csv_meas.local0;
                    meas.local_variance()[0] = csv_meas.var_local0;
                    indices[meas.dimensions()++] = ipar;
                }; break;
                case e_bound_loc1: {
                    meas.local_position()[1] = csv_meas.local1;
                    meas.local_variance()[1] = csv_meas.var_local1;
                    indices[meas.dimensions()++] = ipar;
                }; break;
            }
        }
    }

    meas.time() = csv_meas.time;

    if (acts_to_detray_id) {
        meas.surface_link() = detray::geometry::identifier{
            acts_to_detray_id->at(csv_meas.geometry_id)};
    } else {
        meas.surface_link() =
            detray::geometry::identifier{csv_meas.geometry_id};
    }
    if (det_desc != nullptr) {
        std::size_t dd_idx = geometry_id_to_detector_description_index->at(
            meas.surface_link().value());
        meas.dimensions() = det_desc->dimensions().at(dd_idx);
        meas.set_subspace(det_desc->subspace().at(dd_idx));
    } else {
        meas.set_subspace(indices);
    }
}

}  // namespace traccc::io::csv

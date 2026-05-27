/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_collection.hpp"
#include "traccc/edm/track_state_collection.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

namespace traccc::edm {

/// Types collecting all the collections used to describe (fitted) tracks
template <typename ALGEBRA>
struct track_container {

    struct host {
        /// Constructor using a memory resource
        explicit host(vecmem::memory_resource& mr,
                      measurement_collection::const_view meas = {})
            : tracks{mr}, states{mr}, measurements{meas} {}

        /// The tracks
        track_collection<ALGEBRA>::host tracks;
        /// The track states used by the tracks
        track_state_collection<ALGEBRA>::host states;
        /// The measurements used by the tracks
        measurement_collection::const_view measurements;
    };

    struct buffer {
        /// The tracks
        track_collection<ALGEBRA>::buffer tracks;
        /// The track states used by the tracks
        track_state_collection<ALGEBRA>::buffer states;
        /// The measurements used by the tracks
        measurement_collection::const_view measurements;
    };

    struct data {
        /// Constructor from a host container
        explicit data(host& h)
            : tracks{vecmem::get_data(h.tracks)},
              states{vecmem::get_data(h.states)},
              measurements{h.measurements} {}
        /// The tracks
        track_collection<ALGEBRA>::data tracks;
        /// The track states used by the tracks
        track_state_collection<ALGEBRA>::data states;
        /// The measurements used by the tracks
        measurement_collection::const_view measurements;
    };

    struct const_data {
        /// Constructor from a host container
        explicit const_data(const host& h)
            : tracks{vecmem::get_data(h.tracks)},
              states{vecmem::get_data(h.states)},
              measurements{h.measurements} {}
        /// The tracks
        track_collection<ALGEBRA>::const_data tracks;
        /// The track states used by the tracks
        track_state_collection<ALGEBRA>::const_data states;
        /// The measurements used by the tracks
        measurement_collection::const_view measurements;
    };

    struct view {
        /// Constructor from a buffer
        TRACCC_HOST_DEVICE
        view(const buffer& b)
            : tracks{b.tracks},
              states{b.states},
              measurements{b.measurements} {}
        /// Constructor from a data object
        TRACCC_HOST_DEVICE
        view(const data& d)
            : tracks{d.tracks},
              states{d.states},
              measurements{d.measurements} {}
        /// The tracks
        track_collection<ALGEBRA>::view tracks;
        /// The track states used by the tracks
        track_state_collection<ALGEBRA>::view states;
        /// The measurements used by the tracks
        measurement_collection::const_view measurements;
    };

    struct const_view {
        /// Constructor from a buffer
        TRACCC_HOST_DEVICE
        const_view(const buffer& b)
            : tracks{b.tracks},
              states{b.states},
              measurements{b.measurements} {}
        /// Constructor from a const_data object
        TRACCC_HOST_DEVICE
        const_view(const const_data& d)
            : tracks{d.tracks},
              states{d.states},
              measurements{d.measurements} {}
        /// The tracks
        track_collection<ALGEBRA>::const_view tracks;
        /// The track states used by the tracks
        track_state_collection<ALGEBRA>::const_view states;
        /// The measurements used by the tracks
        measurement_collection::const_view measurements;
    };

    struct device {
        /// Constructor from a view
        TRACCC_HOST_DEVICE
        explicit device(const view& v)
            : tracks{v.tracks},
              states{v.states},
              measurements{v.measurements} {}
        /// The tracks
        track_collection<ALGEBRA>::device tracks;
        /// The track states used by the tracks
        track_state_collection<ALGEBRA>::device states;
        /// The measurements used by the tracks
        measurement_collection::const_device measurements;
    };

    struct const_device {
        /// Constructor from a view
        TRACCC_HOST_DEVICE
        explicit const_device(const const_view& v)
            : tracks{v.tracks},
              states{v.states},
              measurements{v.measurements} {}
        /// The tracks
        track_collection<ALGEBRA>::const_device tracks;
        /// The track states used by the tracks
        track_state_collection<ALGEBRA>::const_device states;
        /// The measurements used by the tracks
        measurement_collection::const_device measurements;
    };

};  // struct track_container

}  // namespace traccc::edm

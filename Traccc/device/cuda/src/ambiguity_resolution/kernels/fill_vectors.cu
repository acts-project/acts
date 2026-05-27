/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "fill_vectors.cuh"

// Project include(s).
#include "traccc/ambiguity_resolution/ambiguity_resolution_config.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>

namespace traccc::cuda::kernels {

__global__ void fill_vectors(const ambiguity_resolution_config cfg,
                             device::fill_vectors_payload payload) {

    const edm::track_container<default_algebra>::const_device track_candidates(
        payload.tracks_view);

    const auto globalIndex = details::global_index1();
    if (globalIndex >= track_candidates.tracks.size()) {
        return;
    }

    const auto track = track_candidates.tracks.at(globalIndex);

    vecmem::jagged_device_vector<measurement_id_type> meas_ids(
        payload.meas_ids_view);
    vecmem::device_vector<measurement_id_type> flat_meas_ids(
        payload.flat_meas_ids_view);
    vecmem::device_vector<traccc::scalar> pvals(payload.pvals_view);
    vecmem::device_vector<unsigned int> n_meas(payload.n_meas_view);
    vecmem::device_vector<int> status(payload.status_view);

    pvals.at(globalIndex) = track.pval();

    if (track.constituent_links().size() < cfg.min_meas_per_track) {
        status.at(globalIndex) = 0;
    } else {
        for (const auto& [type, meas_idx] :
             track_candidates.tracks.constituent_links().at(globalIndex)) {
            assert(type == edm::track_constituent_link::measurement);
            meas_ids.at(globalIndex)
                .push_back(
                    track_candidates.measurements.at(meas_idx).identifier());
            flat_meas_ids.push_back(
                track_candidates.measurements.at(meas_idx).identifier());
        }
        n_meas.at(globalIndex) = track.constituent_links().size();
    }
}
}  // namespace traccc::cuda::kernels

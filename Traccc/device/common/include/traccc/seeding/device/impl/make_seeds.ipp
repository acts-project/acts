#pragma once

#include "traccc/device/global_index.hpp"


namespace traccc::device {

TRACCC_HOST_DEVICE
inline scalar rounded_cot_theta(const scalar zA, const scalar rA,
                                const scalar zB, const scalar rB) {
    static constexpr scalar granularity = 20.f;
    const scalar cot_theta = (zB - zA) / (rB - rA);
    return std::round(cot_theta * granularity) / granularity;
}

TRACCC_HOST_DEVICE
void make_seeds(const global_index_t globalIndex,
                const unsigned int n_images,
                const edm::spacepoint_collection::const_view& spacepoints_view,
                const vecmem::data::vector_view<scalar>& spacepoint_r_view,
                const vecmem::data::vector_view<unsigned int>& nn_outcome_view,
                const vecmem::data::jagged_vector_view<global_index_t>& hits_view,
                vecmem::data::jagged_vector_view<scalar>& cot_theta_matrix_view) {
    if (globalIndex >= n_images) return;
    
    // does not pass the NN selection
    const vecmem::device_vector<unsigned int> nn_outcome(nn_outcome_view);
    if (nn_outcome[globalIndex] == 0) return;
 
    const vecmem::jagged_device_vector<global_index_t> hits(hits_view);
    if (hits[globalIndex].size() < 3) return;

    const edm::spacepoint_collection::const_device spacepoints(spacepoints_view);
    const vecmem::device_vector<scalar> spacepoint_r(spacepoint_r_view);
 
    vecmem::jagged_device_vector<scalar> cot_theta_matrix(cot_theta_matrix_view);

    unsigned int counter = 0u;
    const auto& associated_hits = hits[globalIndex];
    for (unsigned int i(0); i<associated_hits.size(); ++i) {
        const unsigned int i_index = associated_hits[i];
        const scalar ri = spacepoint_r[i_index];
        const scalar zi = spacepoints.z()[i_index];
        for (unsigned int j(0); j<associated_hits.size(); ++j) {
            ++counter;
            if (i == j) continue;
            const unsigned int j_index = associated_hits[j];
            const scalar rj = spacepoint_r[j_index];
            const scalar zj = spacepoints.z()[j_index];
            const scalar cot = rounded_cot_theta(zi, ri, zj, rj);;
            if (cot != cot) continue;
            if (cot == scalar{0} && (zi > 550.f || zi < -550.f)) continue;
            cot_theta_matrix[globalIndex][counter - 1] = cot;
        }
    }
}

} // namespace

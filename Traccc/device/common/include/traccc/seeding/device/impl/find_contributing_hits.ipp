#pragma once

#include "traccc/device/global_index.hpp"


namespace traccc::device {

TRACCC_HOST_DEVICE
void find_contributing_hits(const global_index_t globalIndex,
                            const global_index_t n_cols,
                            const global_index_t n_rows,
                            const unsigned int n_images,
                            const vecmem::data::vector_view<const scalar>& edges_x_view,
                            const vecmem::data::vector_view<const scalar>& edges_y_view,
                            const vecmem::data::vector_view<const scalar>& spacepoint_r_view,
                            const vecmem::data::vector_view<const scalar>& spacepoint_phi_view,
                            const vecmem::data::vector_view<const scalar>& spacepoint_bfield_z_view,
                            const vecmem::data::vector_view<const unsigned long long>& spacepoint_wedge_association_view,
                            const vecmem::data::jagged_vector_view<const global_index_t>& global_index_view,
                            const vecmem::data::vector_view<const unsigned int>& nn_outcome_view,
                            vecmem::data::jagged_vector_view<global_index_t>& hits_view) {
    if (globalIndex >= n_images) return;
    
    // does not pass the NN selection
    const vecmem::device_vector<const unsigned int> nn_outcome(nn_outcome_view);
    if (nn_outcome[globalIndex] == 0) return;
    
    // get coordinates
    const vecmem::jagged_device_vector<const global_index_t> global_index(global_index_view);
    const auto [wedge_idx, element_idx] = get_wedge_element<global_index_t>(globalIndex, global_index);   
 
    const vecmem::device_vector<const scalar> edges_x(edges_x_view);
    const vecmem::device_vector<const scalar> edges_y(edges_y_view);

    const global_index_t peak_global_index = global_index[wedge_idx][element_idx];
    const auto [col, row] = get_col_row(peak_global_index, n_cols, n_rows);
    const scalar min_x = edges_x[col];
    const scalar max_x = edges_x[col + 1];
    const scalar min_y = edges_y[row];
    const scalar max_y = edges_y[row + 1];

   
    vecmem::jagged_device_vector<global_index_t> hits(hits_view);

    const vecmem::device_vector<const unsigned long long> spacepoint_wedge_association(spacepoint_wedge_association_view);
    const vecmem::device_vector<const scalar> spacepoint_r(spacepoint_r_view);
    const vecmem::device_vector<const scalar> spacepoint_phi(spacepoint_phi_view);
    const vecmem::device_vector<const scalar> spacepoint_bfield_z(spacepoint_bfield_z_view);

    for (global_index_t i(0); i<spacepoint_r.capacity(); ++i) {
        // contributes to wedge
        if ((spacepoint_wedge_association[i] & (1ULL << wedge_idx)) == 0ULL) continue;
        scalar min_value_y = transform(min_x, spacepoint_r[i], spacepoint_phi[i], spacepoint_bfield_z[i]);
        scalar max_value_y = transform(max_x, spacepoint_r[i], spacepoint_phi[i], spacepoint_bfield_z[i]);
        if (min_value_y > max_value_y) std::swap(min_value_y, max_value_y);
       
        if (min_value_y > max_y) continue;
        if (max_value_y < min_y) continue;
 
        // add contribution
        hits[globalIndex].push_back(i);
    }
}

} // namespace

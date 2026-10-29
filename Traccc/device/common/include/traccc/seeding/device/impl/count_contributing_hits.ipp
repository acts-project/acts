#pragma once

#include "traccc/device/global_index.hpp"


namespace traccc::device {

template <typename bfield_t,
          typename T,
          unsigned int N_COLS,
          unsigned int N_ROWS>
TRACCC_HOST_DEVICE
void count_contributing_hits(const global_index_t globalIndex,
                             const unsigned int n_images,
                             const bfield_t& bfield,
                             const typename edm::houghplane_collection<T, N_COLS, N_ROWS>::const_view& houghplanes_view,
                             const edm::spacepoint_collection::const_view& spacepoints_view,
                             const vecmem::data::vector_view<scalar>& spacepoint_r_view,
                             const vecmem::data::vector_view<scalar>& spacepoint_phi_view,
                             const vecmem::data::vector_view<unsigned long long>& spacepoint_wedge_association_view,
                             const vecmem::data::jagged_vector_view<global_index_t>& global_index_view,
                             const vecmem::data::vector_view<unsigned int>& nn_outcome_view,
                             vecmem::data::vector_view<unsigned int>& count_vector_view) {
    if (globalIndex >= n_images) return;

    // does not pass the NN selection
    const vecmem::device_vector<unsigned int> nn_outcome(nn_outcome_view);
    if (nn_outcome[globalIndex] == 0) return;
    
    // get coordinates
    const vecmem::jagged_device_vector<global_index_t> global_index(global_index_view);
    
    global_index_t current_size = globalIndex;
    global_index_t wedge_idx = 0;
    global_index_t element_idx = 0;

    for (global_index_t i(0); i<global_index.size(); ++i) {
        if (global_index[i].size() > current_size) {
           wedge_idx = i;
           element_idx = current_size;
           break;
        }
        current_size -= global_index[i].size();
    }
    
    const global_index_t peak_global_index = global_index[wedge_idx][element_idx];
    const typename edm::houghplane_collection<T, N_COLS, N_ROWS>::const_device houghplanes(houghplanes_view);
    const auto [col, row] = get_col_row<N_COLS, N_ROWS>(peak_global_index);
    const scalar min_x = houghplanes.edges_x()[col];
    const scalar max_x = houghplanes.edges_x()[col + 1];
    const scalar min_y = houghplanes.edges_y()[row];
    const scalar max_y = houghplanes.edges_y()[row + 1];

   
    vecmem::device_vector<unsigned int> count_vector(count_vector_view);
    unsigned int n_contributing_hits = 0u;

    const edm::spacepoint_collection::const_device spacepoints(spacepoints_view);
    const vecmem::device_vector<unsigned long long> spacepoint_wedge_association(spacepoint_wedge_association_view);
    const vecmem::device_vector<scalar> spacepoint_r(spacepoint_r_view);
    const vecmem::device_vector<scalar> spacepoint_phi(spacepoint_phi_view);
/*
    printf("We are looking at wedge %u [element %u] from a %ux%u plane, peak is at index %u (i.e. %u, %u) -> entries %hu -> ranges x=%.4f, %.4f y=%.4f, %.4f\n", 
            wedge_idx, element_idx,
            N_COLS, N_ROWS, peak_global_index,
            col, row, 
            houghplanes.plane()[wedge_idx][peak_global_index],
            min_x, max_x,
            min_y, max_y);
*/
    for (unsigned int i(0); i<spacepoints.size(); ++i) {
        // contributes to wedge
        if ((spacepoint_wedge_association[i] & (1ULL << wedge_idx)) == 0ULL) continue;
        const auto covfie_field_at_sp = bfield.at(spacepoints.x()[i], spacepoints.y()[i], spacepoints.z()[i]);

        scalar min_value_y = transform(min_x, spacepoint_r[i], spacepoint_phi[i], covfie_field_at_sp[2]);
        scalar max_value_y = transform(max_x, spacepoint_r[i], spacepoint_phi[i], covfie_field_at_sp[2]);
        
        if (min_value_y > max_y) continue;
        if (max_value_y < min_y) continue;
        
        // add contribution
        ++n_contributing_hits;
    }

    // insert the value
    count_vector[globalIndex] = n_contributing_hits;
}

} // namespace

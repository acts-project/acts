#pragma once

#include "traccc/device/global_index.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"


namespace traccc::device {

template<typename T,
	unsigned int N_COLS,
	unsigned int N_ROWS>
TRACCC_HOST_DEVICE
void ranges(const typename edm::houghplane_collection<T, N_COLS, N_ROWS>::device& container,
	        const scalar magnetic_field_z,
            const scalar r,
            const scalar phi,
            unsigned int* min_ranges,
            unsigned int* max_ranges) {
    for (unsigned int i(0); i<N_COLS; ++i) {
        scalar min_value_phi = transform(container.edges_x()[i], r, phi, magnetic_field_z);
        scalar max_value_phi = transform(container.edges_x()[i + 1], r, phi, magnetic_field_z);
        if (min_value_phi > max_value_phi) std::swap(min_value_phi, max_value_phi);
        
        // out of bounds
        if (min_value_phi < container.edges_y()[0] ||
            max_value_phi > container.edges_y()[N_ROWS]) {
            min_ranges[i] = 0;
            max_ranges[i] = 0;
            return;
        }
        
        unsigned int phi_bin_min = quant(container.edges_y()[0], 
                                         container.edges_y()[N_ROWS],
                                         N_ROWS, 
                                         min_value_phi);
        unsigned int phi_bin_max = quant(container.edges_y()[0], 
                                         container.edges_y()[N_ROWS],
                                         N_ROWS, 
                                         max_value_phi) + 1;
        
        // Clamp bins
        if (phi_bin_max > device::N_ROWS) phi_bin_max = N_ROWS;
        min_ranges[i] = phi_bin_min;
        max_ranges[i] = phi_bin_max;
    }
}


template<typename bfield_t,
    typename T,
    unsigned int N_COLS,
    unsigned int N_ROWS>
TRACCC_HOST_DEVICE
void fill_hough_plane(const global_index_t globalIndex,
                const bfield_t& bfield,
                const edm::spacepoint_collection::const_view::size_type n_spacepoints,
                const edm::spacepoint_collection::const_view& spacepoints_view,
                const vecmem::data::vector_view<unsigned long long>& spacepoint_wedge_association_view,
                const vecmem::data::vector_view<unsigned long long>& spacepoint_detector_id_view,
                const vecmem::data::vector_view<scalar>& spacepoint_r_view,
                const vecmem::data::vector_view<scalar>& spacepoint_phi_view,
                typename edm::houghplane_collection<T, N_COLS, N_ROWS>::view& houghplanes_view) {
    if (globalIndex > N_COLS * N_ROWS) return;

    typename edm::houghplane_collection<T, N_COLS, N_ROWS>::device houghplanes(houghplanes_view);
    const edm::spacepoint_collection::const_device spacepoints(spacepoints_view);

    const vecmem::device_vector<unsigned long long> spacepoint_wedge_association(spacepoint_wedge_association_view);
    const vecmem::device_vector<unsigned long long> spacepoint_detector_id(spacepoint_detector_id_view);
    const vecmem::device_vector<scalar> spacepoint_r(spacepoint_r_view);
    const vecmem::device_vector<scalar> spacepoint_phi(spacepoint_phi_view);

    const auto [x_bin, y_bin] = get_col_row<N_COLS, N_ROWS>(globalIndex);

    // given the bin, take the edges in both dimensions
    const scalar min_x = houghplanes.edges_x()[x_bin];
    const scalar max_x = houghplanes.edges_x()[x_bin + 1];
    const scalar min_y = houghplanes.edges_y()[y_bin];
    const scalar max_y = houghplanes.edges_y()[y_bin + 1];

    std::array<unsigned long long, 25> binContents {};
    for (unsigned int i(0); i<n_spacepoints; ++i) {
      if (spacepoint_wedge_association[i] == 0ULL) continue;

      const auto covfie_field_at_sp = bfield.at(spacepoints.x()[i], spacepoints.y()[i], spacepoints.z()[i]);
      scalar min_value_y = transform(min_x, spacepoint_r[i], spacepoint_phi[i], covfie_field_at_sp[2]);
      scalar max_value_y = transform(max_x, spacepoint_r[i], spacepoint_phi[i], covfie_field_at_sp[2]);

      if (min_value_y > max_y) continue;
      if (max_value_y < min_y) continue;

      for (unsigned int wedge_idx(0); wedge_idx<25; ++wedge_idx) {
          if ((spacepoint_wedge_association[i] & (1ULL << wedge_idx)) == 0ULL) continue;
          binContents[wedge_idx] |= spacepoint_detector_id[i];
      }
    }

    for (unsigned int wedge_idx(0); wedge_idx<25; ++wedge_idx) {
        if (binContents[wedge_idx] == 0ULL) continue;
	    unsigned short value = count_bits(binContents[wedge_idx]);
	    houghplanes.plane()[wedge_idx][globalIndex] = value;
    }
}

} // namespace

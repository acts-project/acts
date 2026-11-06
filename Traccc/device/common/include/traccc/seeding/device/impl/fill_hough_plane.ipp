#pragma once

#include "traccc/device/global_index.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"


namespace traccc::device {

TRACCC_HOST_DEVICE
inline std::pair<unsigned int, unsigned int> range(const global_index_t col_idx,
           const vecmem::device_vector<scalar>& edges_x,
           const vecmem::device_vector<scalar>& edges_y,
	       const scalar magnetic_field_z,
           const scalar r,
           const scalar phi) {
    const unsigned int n_rows = edges_y.size() - 1;

    scalar min_value_phi = transform(edges_x[col_idx], r, phi, magnetic_field_z);
    scalar max_value_phi = transform(edges_x[col_idx + 1], r, phi, magnetic_field_z);
    if (min_value_phi > max_value_phi) std::swap(min_value_phi, max_value_phi);
        
    // out of bounds
    if (min_value_phi < edges_y[0] ||
        max_value_phi > edges_y[n_rows]) {
        return {0u, 0u};
    }
        
    unsigned int phi_bin_min = quant(edges_y[0], 
                                     edges_y[n_rows],
                                     n_rows, 
                                     min_value_phi);
    unsigned int phi_bin_max = quant(edges_y[0], 
                                     edges_y[n_rows],
                                     n_rows, 
                                     max_value_phi) + 1;
        
    // Clamp bins
    if (phi_bin_max > n_rows) phi_bin_max = n_rows;
    return {phi_bin_min, phi_bin_max};  
}

template <typename bfield_t,
          typename T>
TRACCC_HOST_DEVICE
void fill_hough_plane(const global_index_t globalIndex,
		    const bfield_t& bfield,
            const edm::spacepoint_collection::const_view::size_type n_spacepoints,
            const edm::spacepoint_collection::const_view& spacepoints_view,
            const vecmem::data::vector_view<unsigned long long>& spacepoint_wedge_association_view,
            const vecmem::data::vector_view<unsigned long long>& spacepoint_detector_id_view,
		    const vecmem::data::vector_view<scalar>& spacepoint_r_view,
		    const vecmem::data::vector_view<scalar>& spacepoint_phi_view,
            const vecmem::data::vector_view<scalar>& edges_x_view,
            const vecmem::data::vector_view<scalar>& edges_y_view,
            vecmem::data::vector_view<unsigned long long>& bins_view) {
   if (globalIndex >= n_spacepoints) return;
   
   const edm::spacepoint_collection::const_device spacepoints(spacepoints_view);
   const vecmem::device_vector<unsigned long long> spacepoint_wedge_association(spacepoint_wedge_association_view);
   const vecmem::device_vector<unsigned long long> spacepoint_detector_id(spacepoint_detector_id_view);
   const vecmem::device_vector<scalar> spacepoint_r(spacepoint_r_view);
   const vecmem::device_vector<scalar> spacepoint_phi(spacepoint_phi_view);
   vecmem::device_vector<unsigned long long> bins(bins_view);

   // skip if it does not contribute
   if (spacepoint_wedge_association[globalIndex] == 0ULL) return;

   const vecmem::device_vector<scalar> edges_x(edges_x_view);
   const vecmem::device_vector<scalar> edges_y(edges_y_view);
   const global_index_t n_cols = edges_x.size() - 1;
   const global_index_t n_rows = edges_y.size() - 1;

   // get magnetic field
   const auto covfie_field_at_sp = bfield.at(spacepoints.x()[globalIndex], spacepoints.y()[globalIndex], spacepoints.z()[globalIndex]);

   unsigned long long slice = spacepoint_wedge_association[globalIndex];
   const unsigned long long layer_bit = spacepoint_detector_id[globalIndex];

   // Fill the planes
   for (unsigned int x_idx(0); x_idx < n_cols; ++x_idx) {
        const auto [first, last] = range(x_idx, 
                                         edges_x, 
                                         edges_y,
                                         covfie_field_at_sp[2],
                                         spacepoint_r[globalIndex],
                                         spacepoint_phi[globalIndex]);
        if (last == first) continue;
        assert(last > first);
        const unsigned int diff_phi = last - first;

        // loop on wedges
        for (unsigned int wedge_idx(0); wedge_idx<25u; ++wedge_idx) {
            if ((slice & (1ULL << wedge_idx)) == 0ULL) continue;
            const unsigned long long offset = wedge_idx * n_cols * n_rows + x_idx * n_rows;

            if (diff_phi < n_rows / 2) {
                for (unsigned int y_idx(first); y_idx < last; ++y_idx) {
     		        vecmem::device_atomic_ref<unsigned long long> a2(bins[offset + y_idx]);
		            a2.fetch_or(layer_bit);
                }
            } else {
                for (unsigned int y_idx(0); y_idx < first; ++y_idx) {
		            vecmem::device_atomic_ref<unsigned long long> a2(bins[offset + y_idx]);
                    a2.fetch_or(layer_bit);
                }
                for (unsigned int y_idx(last); y_idx < device::N_ROWS; ++y_idx) {
		            vecmem::device_atomic_ref<unsigned long long> a2(bins[offset + y_idx]);
                    a2.fetch_or(layer_bit);
                }
            } // checks on diff_phi
        } // loop on planes
   } // loop on x bins
}

} // namespace

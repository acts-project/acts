#pragma once

#include "traccc/device/global_index.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"


namespace traccc::device {

template <typename T>
TRACCC_HOST_DEVICE
void fill_hough_plane(const global_index_t globalIndex,
            const global_index_t n_wedges,
            const global_index_t n_cols,
            const global_index_t n_rows,
            const edm::spacepoint_collection::const_view::size_type n_spacepoints,
            const vecmem::data::vector_view<const unsigned long long>& spacepoint_wedge_association_view,
            const vecmem::data::vector_view<const unsigned long long>& spacepoint_detector_id_view,
		    const vecmem::data::vector_view<const scalar>& spacepoint_r_view,
		    const vecmem::data::vector_view<const scalar>& spacepoint_phi_view,
            const vecmem::data::vector_view<const scalar>& spacepoint_bfield_z_view,
            const vecmem::data::vector_view<const scalar>& edges_x_view,
            const vecmem::data::vector_view<const scalar>& edges_y_view,
            vecmem::data::vector_view<unsigned long long>& bins_view) {
   if (globalIndex >= n_spacepoints) return;
   
   const vecmem::device_vector<const unsigned long long> spacepoint_wedge_association(spacepoint_wedge_association_view);
   const vecmem::device_vector<const unsigned long long> spacepoint_detector_id(spacepoint_detector_id_view);
   const vecmem::device_vector<const scalar> spacepoint_r(spacepoint_r_view);
   const vecmem::device_vector<const scalar> spacepoint_phi(spacepoint_phi_view);
   const vecmem::device_vector<const scalar> spacepoint_bfield_z(spacepoint_bfield_z_view);
   vecmem::device_vector<unsigned long long> bins(bins_view);

   // skip if it does not contribute
   if (spacepoint_wedge_association[globalIndex] == 0ULL) return;

   const vecmem::device_vector<const scalar> edges_x(edges_x_view);
   const vecmem::device_vector<const scalar> edges_y(edges_y_view);

   unsigned long long slice = spacepoint_wedge_association[globalIndex];
   const unsigned long long layer_bit = spacepoint_detector_id[globalIndex];

   // Fill the planes
   
   // re-use the computation of the transform from the previous bin (right -> left edge)
   scalar left_ht = transform(edges_x[0], spacepoint_r[globalIndex], spacepoint_phi[globalIndex], spacepoint_bfield_z[globalIndex]);
   
   for (unsigned int x_idx(0); x_idx < n_cols; ++x_idx) {
      const scalar right_ht = transform(edges_x[x_idx + 1], spacepoint_r[globalIndex], spacepoint_phi[globalIndex], spacepoint_bfield_z[globalIndex]);
      scalar lo_ht = left_ht;
      scalar hi_ht = right_ht;
      left_ht = right_ht;
      if (lo_ht > hi_ht) std::swap(lo_ht, hi_ht);

      // out of bounds
      if (lo_ht < edges_y[0] ||
          hi_ht > edges_y[n_rows]) continue;

      unsigned int phi_bin_min = quant(edges_y[0],
                                       edges_y[n_rows],
                                       n_rows,
                                       lo_ht);
      unsigned int phi_bin_max = quant(edges_y[0],
                                       edges_y[n_rows],
                                       n_rows,
                                       hi_ht) + 1;

      // Clamp bins
      if (phi_bin_max > n_rows) phi_bin_max = n_rows;
      const unsigned int first = phi_bin_min;
      const unsigned int last = phi_bin_max;
      if (last == first) continue;
      assert(last > first);
      const unsigned int diff_phi = last - first;

        // loop on wedges
        for (unsigned int wedge_idx(0); wedge_idx<n_wedges; ++wedge_idx) {
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
                for (unsigned int y_idx(last); y_idx < n_rows; ++y_idx) {
		            vecmem::device_atomic_ref<unsigned long long> a2(bins[offset + y_idx]);
                    a2.fetch_or(layer_bit);
                }
            } // checks on diff_phi
        } // loop on planes
   } // loop on x bins
}

} // namespace

#pragma once

#include "traccc/device/global_index.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "Acts/Definitions/Units.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
int quant(scalar min, scalar max, unsigned nSteps, scalar val) {
  return static_cast<int>((val - min) / (max - min) * nSteps);
}

TRACCC_HOST_DEVICE
scalar transform(const scalar qopt_value,
                 const scalar r,
                 const scalar phi) {
    // Need to retrieve it from somewhere
    scalar magnetic_field_z = 2.f * Acts::UnitConstants::T;
	    
    const scalar d0 = 0;
    scalar phi_value = asin(0.50f * r * magnetic_field_z * qopt_value - d0 / r) + phi;
    ::remainderf(phi_value, 2.00f * std::numbers::pi_v<scalar>);
    return phi_value;
}

TRACCC_HOST_DEVICE
void ranges(const edm::houghplane_collection<unsigned long long, device::N_COLS, device::N_ROWS>::device& container,
            const scalar r,
            const scalar phi,
            unsigned int* min_ranges,
            unsigned int* max_ranges) {
    for (unsigned int i(0); i<device::N_COLS; ++i) {
        scalar min_value_phi = transform(container.edges_x()[i], r, phi);
        scalar max_value_phi = transform(container.edges_x()[i + 1], r, phi);
        if (min_value_phi > max_value_phi) std::swap(min_value_phi, max_value_phi);
        
        // out of bounds
        if (min_value_phi < container.edges_y()[0] ||
            max_value_phi > container.edges_y()[device::N_ROWS]) {
            min_ranges[i] = 0;
            max_ranges[i] = 0;
            return;
        }
        
        unsigned int phi_bin_min = quant(container.edges_y()[0], 
                                         container.edges_y()[device::N_ROWS],
                                         device::N_ROWS, 
                                         min_value_phi);
        unsigned int phi_bin_max = quant(container.edges_y()[0], 
                                         container.edges_y()[device::N_ROWS],
                                         device::N_ROWS, 
                                         max_value_phi) + 1;
        
        // Clamp bins
        if (phi_bin_max > device::N_ROWS) phi_bin_max = device::N_ROWS;
        min_ranges[i] = phi_bin_min;
        max_ranges[i] = phi_bin_max;
    }
}


TRACCC_HOST_DEVICE
void fill_hough_plane(const global_index_t globalIndex,
                const edm::spacepoint_collection::const_view::size_type n_spacepoints,
                const edm::spacepoint_collection::const_view& spacepoints_view,
                const vecmem::data::vector_view<unsigned long long>& spacepoint_wedge_association_view,
                const vecmem::data::vector_view<unsigned long long>& spacepoint_detector_id_view,
		const vecmem::data::vector_view<scalar>& spacepoint_r_view,
		const vecmem::data::vector_view<scalar>& spacepoint_phi_view,
                edm::houghplane_collection<unsigned long long, device::N_COLS, device::N_ROWS>::view houghplanes_view) {
   if (globalIndex >= n_spacepoints) return;
   
   const edm::spacepoint_collection::const_device spacepoints(spacepoints_view);
   const vecmem::device_vector<unsigned long long> spacepoint_wedge_association(spacepoint_wedge_association_view);
   const vecmem::device_vector<unsigned long long> spacepoint_detector_id(spacepoint_detector_id_view);
   const vecmem::device_vector<scalar> spacepoint_r(spacepoint_r_view);
   const vecmem::device_vector<scalar> spacepoint_phi(spacepoint_phi_view);
   
   edm::houghplane_collection<unsigned long long, device::N_COLS, device::N_ROWS>::device houghplanes(houghplanes_view);

   // skip if it does not contribute
   if (spacepoint_wedge_association[globalIndex] == 0ULL) return;
 
   // arrays for ranges
   unsigned int min_ranges[device::N_COLS] {};
   unsigned int max_ranges[device::N_COLS] {};
   ranges(houghplanes, 
          spacepoint_r[globalIndex], spacepoint_phi[globalIndex], 
	  min_ranges, max_ranges);

   unsigned long long slice = spacepoint_wedge_association[globalIndex];
   const unsigned long long layer_bit = spacepoint_detector_id[globalIndex];

   // Fill the planes
   for (unsigned int x_idx(0); x_idx < device::N_COLS; ++x_idx) {
        const unsigned int first = min_ranges[x_idx];
        const unsigned int last = max_ranges[x_idx];
        if (last == first) continue;
        assert(last > first);
        const unsigned int diff_phi = last - first;

        // loop on wedges
        for (unsigned int wedge_idx(0); wedge_idx<houghplanes.size(); ++wedge_idx) {
            if ((slice & (1ULL << wedge_idx)) == 0ULL) continue;
            const unsigned long long offset = x_idx * device::N_ROWS;

            if (diff_phi < device::N_ROWS / 2) {
                for (unsigned int y_idx(first); y_idx < last; ++y_idx) {
     		   vecmem::device_atomic_ref<unsigned long long> a2(houghplanes.plane()[wedge_idx][offset + y_idx]);
		   a2.fetch_or(layer_bit);
                }
            } else {
                for (unsigned int y_idx(0); y_idx < first; ++y_idx) {
		   vecmem::device_atomic_ref<unsigned long long> a2(houghplanes.plane()[wedge_idx][offset + y_idx]);
                   a2.fetch_or(layer_bit);
                }
                for (unsigned int y_idx(last); y_idx < device::N_ROWS; ++y_idx) {
		   vecmem::device_atomic_ref<unsigned long long> a2(houghplanes.plane()[wedge_idx][offset + y_idx]);
                   a2.fetch_or(layer_bit);
                }
            } // checks on diff_phi
        } // loop on planes
   } // loop on x bins
}
} // namespace

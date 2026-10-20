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
void ranges(const edm::houghplane_collection<device::N_COLS, device::N_ROWS>::device& container,
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


template <detray::concepts::detector detector_t>
TRACCC_HOST_DEVICE
void fill_hough_plane(const global_index_t globalIndex,
		const detray::detector_view_t<detector_t> det_view,
		const edm::measurement_collection::const_view& measurements_view,
		edm::spacepoint_collection::const_view::size_type n_spacepoints,
		const edm::spacepoint_collection::const_view& spacepoints_view,
		const scalar interactionCenter,
		const scalar interactionWindow,
		const scalar etaWidth,
		edm::houghplane_collection<device::N_COLS, device::N_ROWS>::view houghplanes_view) {
   if (globalIndex >= n_spacepoints) return;
   // we need to pass these from outside
   std::array<scalar, 25> wedges = {
  	-2.88f, -2.64f, -2.40f, -2.16f, -1.92f, -1.68f,
    	-1.44f, -1.20f, -0.96f, -0.72f, -0.48f, -0.24f, 0.00f,
    	0.24f, 0.48f, 0.72f, 0.96f, 1.20f, 1.44f,
    	1.68f, 1.92f, 2.16f, 2.40f, 2.64f, 2.88f
   };
   
   const edm::spacepoint_collection::const_device spacepoints(spacepoints_view);
   edm::houghplane_collection<device::N_COLS, device::N_ROWS>::device houghplanes(houghplanes_view);

   scalar z0_min = interactionCenter - interactionWindow;
   scalar z0_max = interactionCenter + interactionWindow;
    
   scalar r = hypotf(spacepoints[globalIndex].x(), spacepoints[globalIndex].y());
   scalar eta_at_min = asinhf((spacepoints[globalIndex].z() - z0_min) / r);
   scalar eta_at_max = asinhf((spacepoints[globalIndex].z() - z0_max) / r);

   // understand what wedge(s) the space points belongs to
   // wedges define eta ranges and they may overlap
   //
   // this implementation has an upper limit on the number of wedges
   assert(houghplanes.size() == wedges.size());
   unsigned long long slice = 0ULL;
    for (unsigned int wedge_idx(0); wedge_idx<houghplanes.size(); ++wedge_idx) {
        if (eta_at_min < wedges[wedge_idx] - etaWidth) continue;
        if (eta_at_max > wedges[wedge_idx] + etaWidth) continue;
        slice |= (1ULL << wedge_idx);
    } // loop on wedges
    if (slice == 0ULL) return;

    // Create the tracking geometry
    typename detray::detector_device_t<detector_t> det(det_view);

    // Retrieve the corresponding measurement index and, from it, the geo id
    // all of this is to retrieve the volume_id and layer_id of the 
    // space point
    const edm::measurement_collection::const_device measurements(measurements_view);
    const auto& meas = measurements.at(spacepoints[globalIndex].measurement_index_1());

    const detray::geometry::identifier& surface_link = meas.surface_link();
    const auto volume_id = meas.surface_link().volume();
    const auto identifier_id = volume_id; // we need to check if we need more ...

    // make hough transformation
    // result is the same for all planes
    scalar phi = atan2f(spacepoints[globalIndex].y(), spacepoints[globalIndex].y());
    ::remainderf(phi, 2.0f * std::numbers::pi_v<scalar>);
    
    // arrays for ranges
    unsigned int min_ranges[device::N_COLS] {};
    unsigned int max_ranges[device::N_COLS] {};
    ranges(houghplanes, r, phi, min_ranges, max_ranges);

    // Fill the planes
    const unsigned long long layer_bit = 1ULL << identifier_id;
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
     		   vecmem::device_atomic_ref<unsigned long long> a2(houghplanes[wedge_idx].plane()[offset + y_idx]);
		   a2.fetch_or(layer_bit);
                }
            } else {
                for (unsigned int y_idx(0); y_idx < first; ++y_idx) {
		   vecmem::device_atomic_ref<unsigned long long> a2(houghplanes[wedge_idx].plane()[offset + y_idx]);
                   a2.fetch_or(layer_bit);
                }
                for (unsigned int y_idx(last); y_idx < device::N_ROWS; ++y_idx) {
		   vecmem::device_atomic_ref<unsigned long long> a2(houghplanes[wedge_idx].plane()[offset + y_idx]);
                   a2.fetch_or(layer_bit);
                }
            } // checks on diff_phi
        } // loop on planes
    } // loop on x bins
}
} // namespace

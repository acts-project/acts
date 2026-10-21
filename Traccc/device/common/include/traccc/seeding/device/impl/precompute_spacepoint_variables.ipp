#pragma once

#include "traccc/device/global_index.hpp"

namespace traccc::device {

template <detray::concepts::detector detector_t>
TRACCC_HOST_DEVICE
void precompute_spacepoint_variables(const global_index_t globalIndex,
		const detray::detector_view_t<detector_t> det_view,
		const edm::measurement_collection::const_view& measurements_view,
		edm::spacepoint_collection::const_view::size_type n_spacepoints,
                const edm::spacepoint_collection::const_view& spacepoints_view,
		const scalar interactionCenter,
                const scalar interactionWindow,
                const scalar etaWidth,
		vecmem::data::vector_view<unsigned long long>& spacepoint_wedge_association_view,
                vecmem::data::vector_view<unsigned long long>& spacepoint_detector_id_view,
		vecmem::data::vector_view<scalar>& spacepoint_r_view,
		vecmem::data::vector_view<scalar>& spacepoint_phi_view) {
   if (globalIndex >= n_spacepoints) return;
   // we need to pass these from outside
   static const std::array<scalar, 25> wedges = {
  	-2.88f, -2.64f, -2.40f, -2.16f, -1.92f, -1.68f,
    	-1.44f, -1.20f, -0.96f, -0.72f, -0.48f, -0.24f, 0.00f,
    	0.24f, 0.48f, 0.72f, 0.96f, 1.20f, 1.44f,
    	1.68f, 1.92f, 2.16f, 2.40f, 2.64f, 2.88f
   };
   
   // input   
   const edm::spacepoint_collection::const_device spacepoints(spacepoints_view);
   // outputs
   vecmem::device_vector<unsigned long long> spacepoint_wedge_association(spacepoint_wedge_association_view);
   vecmem::device_vector<unsigned long long> spacepoint_detector_id(spacepoint_detector_id_view);
   vecmem::device_vector<scalar> spacepoint_r(spacepoint_r_view);
   vecmem::device_vector<scalar> spacepoint_phi(spacepoint_phi_view);

   scalar z0_min = interactionCenter - interactionWindow;
   scalar z0_max = interactionCenter + interactionWindow;
    
   spacepoint_r[globalIndex] = hypotf(spacepoints[globalIndex].x(), spacepoints[globalIndex].y());
   scalar eta_at_min = asinhf((spacepoints[globalIndex].z() - z0_min) / spacepoint_r[globalIndex]);
   scalar eta_at_max = asinhf((spacepoints[globalIndex].z() - z0_max) / spacepoint_r[globalIndex]);

   spacepoint_phi[globalIndex] = atan2f(spacepoints[globalIndex].y(), spacepoints[globalIndex].y());
   ::remainderf(spacepoint_phi[globalIndex], 2.0f * std::numbers::pi_v<scalar>);

   // understand what wedge(s) the space points belongs to
   // wedges define eta ranges and they may overlap
   //
   // this implementation has an upper limit on the number of wedges
   unsigned long long slice = 0ULL;
   for (unsigned int wedge_idx(0); wedge_idx<wedges.size(); ++wedge_idx) {
      if (eta_at_min < wedges[wedge_idx] - etaWidth) continue;
      if (eta_at_max > wedges[wedge_idx] + etaWidth) continue;
      slice |= (1ULL << wedge_idx);
  } // loop on wedges
  if (slice == 0ULL) return; 
  spacepoint_wedge_association[globalIndex] = slice;  

  // Create the tracking geometry
  typename detray::detector_device_t<detector_t> det(det_view);

  // Retrieve the corresponding measurement index and, from it, the geo id
  // all of this is to retrieve the volume_id and layer_id of the 
  // space point
  const edm::measurement_collection::const_device measurements(measurements_view);
  const auto& meas = measurements.at(spacepoints[globalIndex].measurement_index_1());

  const detray::geometry::identifier& surface_link = meas.surface_link();
  const auto volume_id = meas.surface_link().volume();
  const auto identifier_id = volume_id; // we need to check what we need
  spacepoint_detector_id[globalIndex] = identifier_id;
}
} // namespace

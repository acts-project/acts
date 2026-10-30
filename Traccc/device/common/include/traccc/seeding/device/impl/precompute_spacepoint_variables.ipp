#pragma once

#include "traccc/device/global_index.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
unsigned long long getId(const unsigned int index) {
   unsigned long long shift = 0ULL;
   switch (index) {
	   // negative endcap (pixels)
	   case 3: {shift = 0; break;}
       case 5: {shift = 1; break;}
	   case 7: {shift = 2; break;}
	   case 9: {shift = 3; break;}
	   case 11: {shift = 4; break;}
	   case 13: {shift = 5; break;}
	   case 15: {shift = 6; break;}
           // barrel (pixels)
	   case 18: {shift = 7; break;}
	   case 20: {shift = 8; break;}
	   case 22: {shift = 9; break;}
	   case 24: {shift = 10; break;}
 	   // positive endcap (pixels)
	   case 27: {shift = 11; break;}
	   case 29: {shift = 12; break;}
	   case 31: {shift = 13; break;}
	   case 33: {shift = 14; break;}
	   case 35: {shift = 15; break;}
	   case 37: {shift = 16; break;}
	   case 39: {shift = 17; break;}
	   // negative endcap (short strips)
	   case 45: {shift=18; break;}
	   case 47: {shift=19; break;}
	   case 49: {shift=20; break;}
	   case 51: {shift=21; break;}
	   case 53: {shift=22; break;}
	   case 55: {shift=23; break;}
	   // barrel (short strips)
	   case 58: {shift=24; break;}
	   case 60: {shift=25; break;}
	   case 62: {shift=26; break;}
	   case 64: {shift=27; break;}
	   // positive endcap (short strips)
	   case 67: {shift=28; break;}
	   case 69: {shift=29; break;}
	   case 71: {shift=30; break;}
	   case 73: {shift=31; break;}
	   case 75: {shift=32; break;}
	   case 77: {shift=33; break;}
	   // others
	   default:
		   return 0ULL;
   };
   return 1ULL << shift;
}

TRACCC_HOST_DEVICE
void precompute_spacepoint_variables(const global_index_t globalIndex,
		const edm::measurement_collection::const_view& measurements_view,
		const edm::spacepoint_collection::const_view::size_type n_spacepoints,
        const edm::spacepoint_collection::const_view& spacepoints_view,
		const scalar interactionCenter,
        const scalar interactionWindow,
        const scalar etaWidth,
        const unsigned long long module_veto,
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

   spacepoint_phi[globalIndex] = atan2f(spacepoints[globalIndex].y(), spacepoints[globalIndex].x());
   ::remainderf(spacepoint_phi[globalIndex], 2.0f * std::numbers::pi_v<scalar>);

   // Retrieve the corresponding measurement index and, from it, the geo id
   // all of this is to retrieve the volume_id and layer_id of the 
   // space point
   const edm::measurement_collection::const_device measurements(measurements_view);
   const auto& meas = measurements.at(spacepoints[globalIndex].measurement_index_1());

   const auto volume_id = meas.surface_link().volume();
   const auto identifier_id = getId(volume_id);
   spacepoint_detector_id[globalIndex] = identifier_id;

   // veto layers
   if ((identifier_id & module_veto) != 0ULL) return;

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
}
} // namespace

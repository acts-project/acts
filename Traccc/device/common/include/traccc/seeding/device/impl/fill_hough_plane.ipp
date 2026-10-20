#pragma once

#include "traccc/device/global_index.hpp"

namespace traccc::device {

template <detray::concepts::detector detector_t>
TRACCC_HOST_DEVICE
void fill_hough_plane(const global_index_t globalIndex,
		const detray::detector_view_t<detector_t> det_view,
		edm::spacepoint_collection::const_view::size_type n_spacepoints,
		const edm::spacepoint_collection::const_view& spacepoints_view,
		const scalar interactionCenter,
		const scalar interactionWindow,
		edm::houghplane_collection<device::N_COLS, device::N_ROWS>::view houghplanes_view) {
   if (globalIndex >= n_spacepoints) return;
   const edm::spacepoint_collection::const_device spacepoints(spacepoints_view);
   edm::houghplane_collection<device::N_COLS, device::N_ROWS>::device houghplanes(houghplanes_view);

   float z0_min = interactionCenter - interactionWindow;
   float z0_max = interactionCenter + interactionWindow;
    
   float r = hypotf(spacepoints[globalIndex].x(), spacepoints[globalIndex].y());
   float eta_at_min = asinhf((spacepoints[globalIndex].z() - z0_min) / r);
   float eta_at_max = asinhf((spacepoints[globalIndex].z() - z0_max) / r);
}
} // namespace

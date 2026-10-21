#pragma once

#include "traccc/device/global_index.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
void define_plane_edges(const global_index_t globalIndex,
		edm::houghplane_collection<unsigned long long, device::N_COLS, device::N_ROWS>::view houghplanes_view,
		float stepQopt, float minQoptRange,
		float stepPhi, float minPhi) {
	edm::houghplane_collection<unsigned long long, device::N_COLS, device::N_ROWS>::device houghplanes(houghplanes_view);
	if (globalIndex < houghplanes.edges_x().size()) {
		houghplanes.edges_x()[globalIndex] = minQoptRange + stepQopt * globalIndex;
	}
	if (globalIndex < houghplanes.edges_y().size()) {
		houghplanes.edges_y()[globalIndex] = minPhi + stepPhi * globalIndex;
	}
}

} // namespace

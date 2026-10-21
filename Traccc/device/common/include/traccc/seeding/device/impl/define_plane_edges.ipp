#pragma once

#include "traccc/device/global_index.hpp"

namespace traccc::device {

template<typename T,
	unsigned int N_COLS,
	unsigned int N_ROWS>
TRACCC_HOST_DEVICE
void define_plane_edges(const global_index_t globalIndex,
		typename edm::houghplane_collection<T, N_COLS, N_ROWS>::view houghplanes_view,
		float stepQopt, float minQoptRange,
		float stepPhi, float minPhi) {
	typename edm::houghplane_collection<T, N_COLS, N_ROWS>::device houghplanes(houghplanes_view);
	if (globalIndex < houghplanes.edges_x().size()) {
		houghplanes.edges_x()[globalIndex] = minQoptRange + stepQopt * globalIndex;
	}
	if (globalIndex < houghplanes.edges_y().size()) {
		houghplanes.edges_y()[globalIndex] = minPhi + stepPhi * globalIndex;
	}
}

} // namespace

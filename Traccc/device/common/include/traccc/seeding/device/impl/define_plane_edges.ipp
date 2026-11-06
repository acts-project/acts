#pragma once

#include "traccc/device/global_index.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
void define_plane_edges(const global_index_t globalIndex,
        vecmem::data::vector_view<scalar>& edges_x_view,
        vecmem::data::vector_view<scalar>& edges_y_view,
		const scalar stepQopt, 
        const scalar minQoptRange,
		const scalar stepPhi, 
        const scalar minPhi) {
    vecmem::device_vector<scalar> edges_x(edges_x_view);
    vecmem::device_vector<scalar> edges_y(edges_y_view);
	
    if (globalIndex < edges_x.capacity()) edges_x[globalIndex] = minQoptRange + stepQopt * globalIndex;
    if (globalIndex < edges_y.capacity()) edges_y[globalIndex] = minPhi + stepPhi * globalIndex;
}

} // namespace

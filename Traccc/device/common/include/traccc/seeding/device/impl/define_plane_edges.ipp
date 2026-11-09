#pragma once

#include "traccc/device/global_index.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
void define_plane_edges(const global_index_t globalIndex,
        const global_index_t n_cols,
        const global_index_t n_rows,
		const scalar stepQopt, 
        const scalar minQoptRange,
		const scalar stepPhi, 
        const scalar minPhi,
        vecmem::data::vector_view<scalar>& edges_x_view,
        vecmem::data::vector_view<scalar>& edges_y_view) {
    vecmem::device_vector<scalar> edges_x(edges_x_view);
    vecmem::device_vector<scalar> edges_y(edges_y_view);

    if (globalIndex <= n_cols) edges_x[globalIndex] = minQoptRange + stepQopt * globalIndex;
    if (globalIndex <= n_rows) edges_y[globalIndex] = minPhi + stepPhi * globalIndex;
}

} // namespace

#pragma once

#include "traccc/device/global_index.hpp"

namespace traccc::device {

template <typename T,
          unsigned int N_COLS,
          unsigned int N_ROWS>
TRACCC_HOST_DEVICE
void evaluate_nn(const global_index_t globalIndex,
        const unsigned int n_images,
        const typename edm::houghplane_collection<T, N_ROWS, N_COLS>::view& images_view,
        const vecmem::data::vector_view<scalar>& maxima_x_view,
        vecmem::data::vector_view<unsigned int>& nn_outcome_view) {
    if (globalIndex >= n_images) return;
    
    const typename edm::houghplane_collection<T, N_ROWS, N_COLS>::device images(images_view);
    const vecmem::device_vector<scalar> maxima_x(maxima_x_view);
    vecmem::device_vector<unsigned int> nn_outcome(nn_outcome_view);

    // For now we pass everything
    nn_outcome[globalIndex] = 1u;
}
} // namespace

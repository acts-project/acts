/**
 * TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <cuda_runtime.h>
#include <gtest/gtest.h>

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);

    int r = RUN_ALL_TESTS();

    cudaError_t cErr = cudaDeviceReset();

    if (cErr == cudaSuccess || cErr == cudaErrorNoDevice ||
        cErr == cudaErrorInsufficientDriver) {
        return r;
    } else {
        return static_cast<int>(cErr);
    }
}

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/examples/hough_transform_seeding_cuda.hpp"
#include "traccc/examples/cuda/full_chain_hough_transform_algorithm.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/host_memory_resource.hpp>

int main(int argc, char* argv[]) {
  // We use a pinned CUDA host memory resource to allocate memory for the
  // cell inputs in order to speed up the memory copies.
  vecmem::cuda::host_memory_resource pinned_host_mr;

  // Execute the throughput test.
  return traccc::hough_transform_seeding_cuda<traccc::cuda::full_chain_hough_transform_algorithm>(
      "Multi-threaded CUDA GPU hough transform seeding test", argc, argv, &pinned_host_mr);
}

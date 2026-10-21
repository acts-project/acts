// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/examples/throughput_mt.hpp"

#include "traccc/examples/alpaka/full_chain_algorithm.hpp"

int main(int argc, char* argv[]) {
  // Execute the throughput test.
  return traccc::throughput_mt<traccc::alpaka::full_chain_algorithm>(
      "Multi-threaded Alpaka GPU throughput tests", argc, argv);
}

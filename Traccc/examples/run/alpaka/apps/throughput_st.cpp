/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/examples/throughput_st.hpp"

#include "traccc/examples/alpaka/full_chain_algorithm.hpp"

int main(int argc, char* argv[]) {

    // Execute the throughput test.
    return traccc::throughput_st<traccc::alpaka::full_chain_algorithm>(
        "Single-threaded Alpaka GPU throughput tests", argc, argv);
}

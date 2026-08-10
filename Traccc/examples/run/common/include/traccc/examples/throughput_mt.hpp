/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <string_view>

namespace traccc {

/// Helper function running a multi-threaded throughput test
///
/// @tparam FULL_CHAIN_ALG The type of the full chain algorithm to use
///
/// @param description A short description of the application
/// @param argc The count of command line arguments (from @c main(...))
/// @param argv The command line arguments (from @c main(...))
/// @param host_mr An optional pointer to a (pinned) host memory resource
///
/// @return The value to be returned from @c main(...)
///
template <typename FULL_CHAIN_ALG>
int throughput_mt(std::string_view description, int argc, char* argv[],
                  vecmem::memory_resource* host_mr = nullptr);

}  // namespace traccc

// Local include(s).
#include "traccc/examples/impl/throughput_mt.ipp"

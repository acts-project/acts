/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

namespace traccc {

// Simple struct for combining multiple memory resources
struct memory_resource {

    // device or shared memory resource
    vecmem::memory_resource& main;

    // optional host accesible memory resource
    vecmem::memory_resource* host = nullptr;
};

}  // namespace traccc

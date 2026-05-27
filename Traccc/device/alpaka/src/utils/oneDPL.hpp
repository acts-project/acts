/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Mark this as a "system header". To suppress all warnings from oneDPL.
// This is needed because at the time of writing we cannot provide oneDPL with
// "-isystem" to the oneAPI compiler.
#pragma clang system_header

// oneDPL include(s).
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>

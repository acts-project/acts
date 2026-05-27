/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::cuda::details {

/// Get current CUDA device number.
///
/// This function wraps the cudaGetDevice function in a way that returns the
/// device number rather than use a reference argument to write to.
///
/// Note that calling the function on a machine with no CUDA device does not
/// result in an error, the function just returns -1 in that case.
///
int get_device();

}  // namespace traccc::cuda::details

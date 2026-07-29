/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc {

/// Enumeration of strategies to await for completion of asynchronous
/// operations.
enum class await_strategy {
  sync_event,  ///< Synchronous waiting for an event to complete
  callback     ///< Suspension with a callback function
};

}  // namespace traccc

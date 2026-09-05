/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// VecMem include(s).
#include <vecmem/utils/abstract_event.hpp>

namespace traccc::device {

/// Abstract base class for synchronizing asynchronous operations on a device
///
/// This class provides an interface for synchronizing asynchronous operations
/// on a device. Derived classes must implement the await method to synchronize
/// an event related to asynchronous operations.
class abstract_awaitable {
 public:
  virtual ~abstract_awaitable() = default;

  /// Synchronize an event related to asynchronous operations
  /// @param event The event to synchronize
  ///
  virtual void await(vecmem::abstract_event& event) const = 0;
};  // class abstract_awaitable

}  // namespace traccc::device

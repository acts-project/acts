// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <functional>
#include <memory>

namespace Acts::detail {

/// A minimal fork-join worker pool for data-parallel loops inside Core.
///
/// The pool owns `numThreads - 1` persistent worker threads; the thread
/// calling @ref parallelFor participates as the remaining worker. A pool
/// with `numThreads <= 1` owns no threads and runs every loop serially on
/// the calling thread.
///
/// Calls to @ref parallelFor from different threads are serialized. A call
/// issued from inside a loop body of the same pool runs serially instead of
/// dead-locking.
class WorkerPool {
 public:
  /// Create a pool with the given total concurrency
  /// @param numThreads Total number of threads including the caller
  explicit WorkerPool(std::size_t numThreads);

  WorkerPool(const WorkerPool&) = delete;
  WorkerPool& operator=(const WorkerPool&) = delete;
  WorkerPool(WorkerPool&&) = delete;
  WorkerPool& operator=(WorkerPool&&) = delete;

  /// Stops and joins the worker threads
  ~WorkerPool();

  /// Total concurrency including the calling thread
  /// @return the number of threads used by @ref parallelFor
  std::size_t size() const;

  /// Run @p body for every index in `[0, n)` and block until all are done.
  ///
  /// Indices are handed out dynamically in chunks of @p grainSize, so the
  /// execution order across indices is unspecified. The first exception
  /// thrown by any invocation is rethrown on the calling thread after all
  /// threads have stopped working on the loop.
  ///
  /// @param n Number of loop iterations
  /// @param body Callable invoked with each index
  /// @param grainSize Number of consecutive indices claimed at once,
  ///        0 picks a default based on @p n and the pool size
  void parallelFor(std::size_t n, const std::function<void(std::size_t)>& body,
                   std::size_t grainSize = 0);

 private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

}  // namespace Acts::detail

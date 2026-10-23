// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/detail/WorkerPool.hpp"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace Acts::detail {

namespace {
/// The pool whose loop body the current thread is executing, if any
thread_local const void* tlActivePool = nullptr;
}  // namespace

struct WorkerPool::Impl {
  std::vector<std::thread> threads;

  /// Serializes concurrent parallelFor callers
  std::mutex callMutex;

  std::mutex mutex;
  std::condition_variable wake;
  std::condition_variable done;
  std::size_t generation = 0;
  std::size_t busy = 0;
  bool stop = false;

  // The current loop, written under `mutex` before the generation bump
  const std::function<void(std::size_t)>* body = nullptr;
  std::size_t n = 0;
  std::size_t grainSize = 1;
  std::atomic<std::size_t> next{0};
  std::atomic<bool> failed{false};
  std::exception_ptr error;

  void work() {
    tlActivePool = this;
    while (!failed.load(std::memory_order_relaxed)) {
      std::size_t begin = next.fetch_add(grainSize, std::memory_order_relaxed);
      if (begin >= n) {
        break;
      }
      std::size_t end = std::min(n, begin + grainSize);
      try {
        for (std::size_t i = begin; i < end; ++i) {
          (*body)(i);
        }
      } catch (...) {
        std::lock_guard lock(mutex);
        if (!error) {
          error = std::current_exception();
        }
        failed.store(true, std::memory_order_relaxed);
      }
    }
    tlActivePool = nullptr;
  }

  void workerLoop() {
    std::size_t seen = 0;
    while (true) {
      {
        std::unique_lock lock(mutex);
        wake.wait(lock, [&] { return stop || generation != seen; });
        if (stop) {
          return;
        }
        seen = generation;
      }
      work();
      std::lock_guard lock(mutex);
      if (--busy == 0) {
        done.notify_one();
      }
    }
  }
};

WorkerPool::WorkerPool(std::size_t numThreads)
    : m_impl(std::make_unique<Impl>()) {
  for (std::size_t i = 1; i < numThreads; ++i) {
    m_impl->threads.emplace_back([impl = m_impl.get()] { impl->workerLoop(); });
  }
}

WorkerPool::~WorkerPool() {
  {
    std::lock_guard lock(m_impl->mutex);
    m_impl->stop = true;
  }
  m_impl->wake.notify_all();
  for (auto& thread : m_impl->threads) {
    thread.join();
  }
}

std::size_t WorkerPool::size() const {
  return m_impl->threads.size() + 1;
}

void WorkerPool::parallelFor(std::size_t n,
                             const std::function<void(std::size_t)>& body,
                             std::size_t grainSize) {
  // Serial execution: no workers, trivial loop, or nested call from a body
  if (m_impl->threads.empty() || n <= 1 || tlActivePool == m_impl.get()) {
    for (std::size_t i = 0; i < n; ++i) {
      body(i);
    }
    return;
  }

  std::lock_guard call(m_impl->callMutex);
  Impl& impl = *m_impl;
  {
    std::lock_guard lock(impl.mutex);
    impl.body = &body;
    impl.n = n;
    impl.grainSize =
        grainSize > 0 ? grainSize : std::max<std::size_t>(1, n / (8 * size()));
    impl.next.store(0, std::memory_order_relaxed);
    impl.failed.store(false, std::memory_order_relaxed);
    impl.error = nullptr;
    impl.busy = impl.threads.size();
    ++impl.generation;
  }
  impl.wake.notify_all();

  impl.work();

  std::exception_ptr error;
  {
    std::unique_lock lock(impl.mutex);
    impl.done.wait(lock, [&] { return impl.busy == 0; });
    impl.body = nullptr;
    error = std::exchange(impl.error, nullptr);
  }
  if (error) {
    std::rethrow_exception(error);
  }
}

}  // namespace Acts::detail

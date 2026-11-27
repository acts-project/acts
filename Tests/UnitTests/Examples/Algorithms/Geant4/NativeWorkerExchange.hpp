// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <condition_variable>
#include <cstddef>
#include <exception>
#include <future>
#include <list>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ActsTests::NativeWorkers {

struct Request {
  int eventId;
  long seed;
  double x;
  bool fail = false;
};

struct Result {
  int eventId;
  double random;
  /// Hit positions (x, z) in millimetres.
  std::vector<std::pair<double, double>> hits;
};

struct Job {
  explicit Job(Request input) : request(input) {}

  const Request request;

 private:
  friend class Exchange;
  std::promise<Result> m_response;
};

/// Bounded exchange between ACTS callers and Geant4 workers in the test.
/// Taking a job retains its slot until the worker completes or fails it.
class Exchange {
 public:
  explicit Exchange(std::size_t capacity) : m_capacity(capacity) {
    if (capacity == 0) {
      throw std::invalid_argument("Exchange capacity must be positive");
    }
  }

  std::future<Result> submit(Request request) {
    std::unique_lock lock(m_mutex);
    m_changed.wait(lock, [&] {
      return m_closed || m_ready.size() + m_running.size() < m_capacity;
    });
    if (m_cancelled) {
      std::rethrow_exception(m_cancelled);
    }
    if (m_closed) {
      throw std::runtime_error("Worker exchange is closed");
    }
    auto job = std::make_shared<Job>(request);
    auto result = job->m_response.get_future();
    m_ready.push_back(std::move(job));
    m_peak = std::max(m_peak, m_ready.size() + m_running.size());
    m_changed.notify_all();
    return result;
  }

  std::shared_ptr<Job> take() {
    std::unique_lock lock(m_mutex);
    m_changed.wait(lock, [&] { return m_closed || !m_ready.empty(); });
    if (m_ready.empty()) {
      return nullptr;
    }
    m_running.splice(m_running.end(), m_ready, m_ready.begin());
    return m_running.back();
  }

  void complete(const std::shared_ptr<Job>& job, Result result) {
    finish(job, [&](auto& response) { response.set_value(std::move(result)); });
  }

  void fail(const std::shared_ptr<Job>& job, std::exception_ptr error) {
    if (!error) {
      throw std::invalid_argument("A failed job requires an exception");
    }
    finish(job, [&](auto& response) { response.set_exception(error); });
  }

  /// Finish accepting jobs. Workers can still drain all admitted requests.
  void close() {
    std::lock_guard lock(m_mutex);
    m_closed = true;
    m_changed.notify_all();
  }

  /// Fail every unresolved request, including those held by workers.
  void cancel(const std::exception_ptr& error) {
    if (!error) {
      throw std::invalid_argument("Cancellation requires an exception");
    }
    std::lock_guard lock(m_mutex);
    if (m_cancelled) {
      return;
    }
    m_cancelled = error;
    m_closed = true;
    m_running.splice(m_running.end(), m_ready);
    for (auto& job : m_running) {
      job->m_response.set_exception(error);
    }
    m_running.clear();
    m_changed.notify_all();
  }

  std::size_t maxOutstanding() const {
    std::lock_guard lock(m_mutex);
    return m_peak;
  }

 private:
  template <typename Resolve>
  void finish(const std::shared_ptr<Job>& job, Resolve resolve) {
    std::lock_guard lock(m_mutex);
    auto it = std::ranges::find(m_running, job);
    if (it == m_running.end()) {
      return;  // Cancellation may have resolved the worker's request.
    }
    resolve((*it)->m_response);
    m_running.erase(it);
    m_changed.notify_all();
  }

  const std::size_t m_capacity;
  mutable std::mutex m_mutex;
  std::condition_variable m_changed;
  std::list<std::shared_ptr<Job>> m_ready;
  std::list<std::shared_ptr<Job>> m_running;
  std::size_t m_peak = 0;
  bool m_closed = false;
  std::exception_ptr m_cancelled;
};

}  // namespace ActsTests::NativeWorkers

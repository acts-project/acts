// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/detail/ScopeExit.hpp"
#include "ActsExamples/EventData/SimParticleDecay.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"

#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <exception>
#include <future>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <utility>

namespace ActsExamples::Geant4 {

struct NativeEventRequest {
  NativeEventRequest(const SimParticleContainer& particles,
                     const Acts::GeometryContext& geometry,
                     const SimParticleDecayGraph* decayGraph,
                     std::uint64_t eventSeed, int eventId)
      : store(geometry, particles, decayGraph), seed(eventSeed), id(eventId) {}
  EventStore store;
  const std::uint64_t seed;
  const int id;
  std::promise<EventResult> response;
};

class NativeEventQueue {
 public:
  explicit NativeEventQueue(std::size_t capacity) : m_capacity(capacity) {}

  EventResult submit(const SimParticleContainer& particles,
                     const Acts::GeometryContext& geometry,
                     const SimParticleDecayGraph* decayGraph,
                     std::uint64_t seed, int eventNumber) {
    std::future<EventResult> result;
    {
      std::unique_lock lock(m_mutex);
      m_changed.wait(
          lock, [&] { return m_closed || !m_open || m_inFlight < m_capacity; });
      if (m_error) {
        std::rethrow_exception(m_error);
      }
      if (!m_open || m_closed) {
        throw std::logic_error("Geant4 native run is not accepting events");
      }
      auto request = std::make_shared<NativeEventRequest>(
          particles, geometry, decayGraph, seed, eventNumber);
      result = request->response.get_future();
      m_ready.push_back(std::move(request));
      ++m_inFlight;
    }
    // Keep at most one admitted event per Geant4 worker, including completed
    // results whose submitting ACTS task has not yet received them.
    Acts::detail::ScopeExit release([&] {
      std::lock_guard lock(m_mutex);
      --m_inFlight;
      m_changed.notify_all();
    });
    m_changed.notify_all();
    return result.get();
  }

  std::shared_ptr<NativeEventRequest> take() {
    std::unique_lock lock(m_mutex);
    // No input is admitted until every worker has completed native run
    // initialization and requested its first event.
    if (!m_open && !m_closed && ++m_readyWorkers == m_capacity) {
      m_open = true;
      m_changed.notify_all();
    }
    m_changed.wait(lock, [&] { return m_closed || !m_ready.empty(); });
    if (m_ready.empty()) {
      return nullptr;
    }
    auto request = std::move(m_ready.front());
    m_ready.pop_front();
    return request;
  }

  void waitUntilOpen() {
    std::unique_lock lock(m_mutex);
    m_changed.wait(lock, [&] { return m_open || m_closed; });
    if (m_error) {
      std::rethrow_exception(m_error);
    }
    if (m_closed) {
      throw std::logic_error("Geant4 native run closed during initialization");
    }
  }

  void close() {
    std::lock_guard lock(m_mutex);
    m_closed = true;
    m_changed.notify_all();
  }

  /// Fail unstarted requests immediately. A worker finishes an active request
  /// itself, so its synchronous caller's event data always outlives transport.
  void fail(std::exception_ptr error) {
    std::lock_guard lock(m_mutex);
    if (!m_error) {
      m_error = std::move(error);
    }
    m_closed = true;
    for (auto& request : m_ready) {
      request->response.set_exception(m_error);
    }
    m_ready.clear();
    m_changed.notify_all();
  }

  void checkError() {
    std::lock_guard lock(m_mutex);
    if (m_error) {
      std::rethrow_exception(m_error);
    }
  }

  void checkFinished() {
    std::lock_guard lock(m_mutex);
    if (m_error) {
      std::rethrow_exception(m_error);
    }
    if (!m_closed) {
      throw std::runtime_error(
          "Geant4 native run ended before the event queue closed");
    }
  }

 private:
  const std::size_t m_capacity;
  std::mutex m_mutex;
  std::condition_variable m_changed;
  std::deque<std::shared_ptr<NativeEventRequest>> m_ready;
  std::size_t m_inFlight = 0;
  std::size_t m_readyWorkers = 0;
  bool m_open = false;
  bool m_closed = false;
  std::exception_ptr m_error;
};

}  // namespace ActsExamples::Geant4

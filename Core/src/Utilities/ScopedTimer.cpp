// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/ScopedTimer.hpp"

#include <chrono>
#include <cmath>

namespace Acts {

ScopedTimer::ScopedTimer(const std::string& name, const Logger& logger,
                         Logging::Level lvl)
    : m_name(name), m_lvl(lvl), m_logger(&logger) {
  m_start = clock_type::now();
}

ScopedTimer::~ScopedTimer() {
  auto end = clock_type::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end - m_start);
  if (m_logger->doPrint(m_lvl)) {
    std::ostringstream oss;
    oss << m_name << " took " << (duration.count() * 1e-3) << " ms";
    m_logger->log(m_lvl, oss.str());
  }
}

AveragingScopedTimer::AveragingScopedTimer(const std::string& name,
                                           const Logger& logger,
                                           Logging::Level lvl)
    : m_name(name), m_lvl(lvl), m_logger(&logger) {}

void AveragingScopedTimer::addSample(std::chrono::nanoseconds duration) {
  const double ns = static_cast<double>(duration.count());
  m_sumDuration.fetch_add(ns, std::memory_order_relaxed);
  m_sumDurationSquared.fetch_add(ns * ns, std::memory_order_relaxed);
  m_nSamples.fetch_add(1, std::memory_order_relaxed);
}

AveragingScopedTimer::Sample::Sample(AveragingScopedTimer& parent)
    : m_parent(&parent), m_start(clock_type::now()) {}

AveragingScopedTimer::Sample::Sample(Sample&& other) noexcept
    : m_parent(other.m_parent), m_start(other.m_start) {
  // Disable the moved-from sample
  other.m_parent = nullptr;
}

AveragingScopedTimer::Sample AveragingScopedTimer::sample() {
  return Sample{*this};
}

AveragingScopedTimer::Sample::~Sample() {
  // parent nullptr means the sample was moved
  if (m_parent != nullptr) {
    auto end = clock_type::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - m_start);
    m_parent->addSample(duration);
  }
}

AveragingScopedTimer::~AveragingScopedTimer() {
  if (m_logger->doPrint(m_lvl)) {
    const double sumDuration =
        m_sumDuration.load(std::memory_order_relaxed);
    const double sumDurationSquared =
        m_sumDurationSquared.load(std::memory_order_relaxed);
    const std::size_t nSamples = m_nSamples.load(std::memory_order_relaxed);
    std::ostringstream oss;
    if (nSamples > 0) {
      double mean = sumDuration / nSamples;
      double stddev = std::sqrt(sumDurationSquared / nSamples - mean * mean);
      oss << m_name << " took " << (sumDuration * 1e-6) << " ms total, "
          << (mean * 1e-3) << " us +- " << (stddev * 1e-3)
          << " us per sample (#" << nSamples << ")";
    } else {
      oss << m_name << " took " << (sumDuration * 1e-6)
          << " ms total (no samples)";
    }
    m_logger->log(m_lvl, oss.str());
  }
}

}  // namespace Acts

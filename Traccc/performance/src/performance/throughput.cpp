// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/performance/throughput.hpp"

// System include(s).
#include <chrono>
#include <iomanip>
#include <iostream>

namespace traccc::performance {

throughput::throughput(std::size_t events, const timing_info& ti,
                       std::string_view timer_name)
    : m_timer_name(timer_name) {
  auto totalTime = ti.get_time(timer_name);
  m_perEvent = (totalTime / static_cast<double>(events)).count() * 1e-6;
  m_perSecond = 1000. / m_perEvent;
}

std::ostream& operator<<(std::ostream& out, const throughput& thr) {
  out << std::setw(30) << std::right << thr.m_timer_name << "  "
      << thr.m_perEvent << " ms/event, " << thr.m_perSecond << " events/s";
  return out;
}

}  // namespace traccc::performance

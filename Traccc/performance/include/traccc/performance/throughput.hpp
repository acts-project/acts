// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/performance/timing_info.hpp"

// System include(s).
#include <iosfwd>
#include <string>
#include <string_view>

namespace traccc::performance {

/// Convenience type for printing throughput information
struct throughput {
  /// Constructor with a "timer name" and event processing time
  throughput(std::size_t events, const timing_info& ti,
             std::string_view timer_name);

  /// The timer name
  std::string m_timer_name;
  /// The milliseconds per event value
  double m_perEvent;
  /// The events per second value
  double m_perSecond;

};  // class throughput

/// Printout operator
std::ostream& operator<<(std::ostream& out, const throughput& thr);

}  // namespace traccc::performance

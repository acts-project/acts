// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/options/details/interface.hpp"

// System include(s).
#include <cstddef>

namespace traccc::opts {

/// Option(s) for multi-threaded code execution
class threading : public interface {
 public:
  /// @name Options
  /// @{

  enum class await_strategy {
    sync_event,  ///< Synchronous waiting for an event to complete
    callback     ///< Suspension with a callback function
  };

  /// The strategy to use for awaiting for asynchronous operations to complete
  await_strategy await_mode = await_strategy::sync_event;

  /// The number of threads to use for the data processing
  std::size_t threads = 1;

  /// The number of events that can be processed concurrently
  std::size_t concurrent_slots = 1;

  /// @}

  /// Constructor
  threading();

  /// Read/process the command line options
  ///
  /// @param vm The command line options to interpret/read
  ///
  void read(const boost::program_options::variables_map& vm) override;

  std::unique_ptr<configuration_printable> as_printable() const override;
};  // struct threading

}  // namespace traccc::opts

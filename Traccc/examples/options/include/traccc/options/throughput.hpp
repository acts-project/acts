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
#include <string>

namespace traccc::opts {

/// Command line options used in the throughput tests
class throughput : public interface {
 public:
  /// @name Options
  /// @{

  /// "Reconstruction stage" to run
  enum class stage {
    seeding,  ///< Run until the end of seeding
    full      ///< Run the full chain of reconstruction
  };
  /// The reconstruction stage to run
  stage reco_stage = stage::full;

  /// The number of events to process during the job
  std::size_t processed_events = 100;
  /// The number of events to run "cold", i.e. run without accounting for
  /// them in the performance measurements
  std::size_t cold_run_events = 10;

  /// Enable or disable the randomization of event processing
  bool deterministic_event_order = true;
  /// Set the random event processing seed
  unsigned int random_seed = 0;

  /// Output log file
  std::string log_file;

  /// @}

  /// Constructor
  throughput();

  std::unique_ptr<configuration_printable> as_printable() const override;

  /// Read/process the command line options
  ///
  /// @param vm The command line options to interpret/read
  ///
  void read(const boost::program_options::variables_map& vm) override;

};  // class throughput

}  // namespace traccc::opts

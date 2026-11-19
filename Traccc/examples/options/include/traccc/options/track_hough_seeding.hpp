// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/options/details/config_provider.hpp"
#include "traccc/options/details/interface.hpp"
#include "traccc/seeding/detail/hough_seeding_config.hpp"

// System include(s).
#include <cstddef>

namespace traccc::opts {

/// Option(s) for multi-threaded code execution
class track_hough_seeding : public interface,
                            public config_provider<hough_seeding_config> {
 public:
  /// @name Options
  /// @{
  bool useHough = false;

  /// @}
  // algorithm config
  hough_seeding_config hough_config;
  explicit operator hough_seeding_config() const override;

  /// Constructor
  track_hough_seeding();

  /// Read/process the command line options
  ///
  /// @param vm The command line options to interpret/read
  ///
  void read(const boost::program_options::variables_map& vm) override;

  std::unique_ptr<configuration_printable> as_printable() const override;
};

}  // namespace traccc::opts

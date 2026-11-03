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
#include "traccc/utils/logging.hpp"

// Boost include(s).
#include <boost/program_options.hpp>

// System include(s).
#include <functional>
#include <string_view>
#include <vector>

namespace traccc::opts {

/// Top-level propgram options for an executable
class program_options {
 public:
  /// Constructor
  program_options(
      std::string_view description,
      const std::vector<std::reference_wrapper<interface> >& options, int argc,
      char* argv[],
      std::unique_ptr<const traccc::Logger> ilogger =
          traccc::getDummyLogger().clone());

 private:
  /// Description of all program options
  boost::program_options::options_description m_desc;

};  // class program_options

}  // namespace traccc::opts

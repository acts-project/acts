// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/utils/logging.hpp"

// Boost
#include "detray/options/boost_program_options.hpp"

// System include(s).
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace detray::options {

// Forward declare the options handling for different configuration types
template <typename T>
void add_options(boost::program_options::options_description& /*desc*/,
                 const T& /*opt*/);

template <typename T>
void configure_options(const boost::program_options::variables_map& /*map*/,
                       T& /*opt*/);

/// Parse commandline options and add them to detray configuration types
template <typename... CONFIGS>
auto parse_options(boost::program_options::options_description& desc, int argc,
                   char* argv[], CONFIGS&... cfgs) {
  static_assert(sizeof...(CONFIGS) > 0, "No commandline options configured");

  if (!argv) {
    throw std::invalid_argument("Invalid command line arguments passed");
  }

  desc.add_options()("help", "Produce help message");

  // Add options according to the configurations that were passed
  (add_options(desc, cfgs), ...);

  // Parse options
  boost::program_options::variables_map vm;
  try {
    boost::program_options::store(
        parse_command_line(
            argc, argv, desc,
            boost::program_options::command_line_style::unix_style ^
                boost::program_options::command_line_style::allow_short),
        vm);

    boost::program_options::notify(vm);
  } catch (const std::exception& ex) {
    // Print help message in case of error
    DETRAY_FATAL_HOST(ex.what() << "\n" << desc);
    std::exit(EXIT_FAILURE);
  }

  // Print help message when requested
  if (vm.count("help") != 0u) {
    std::clog << desc << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  // Add the options to the configurations
  (configure_options(vm, cfgs), ...);
  // Make sure everything is configured correctly
  (print_options(cfgs), ...);

  return vm;
}

/// Parse commandline options and add them to detray configuration types
template <typename... CONFIGS>
auto parse_options(const std::string& description, int argc, char* argv[],
                   CONFIGS&... cfgs) {
  // Options description
  boost::program_options::options_description desc(description);

  // Run options parsing
  return parse_options(desc, argc, argv, cfgs...);
}

}  // namespace detray::options

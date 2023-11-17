// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "CommandLineArguments.h"

#include "Acts/Plugins/Sycl/Utilities/ListPlatforms.hpp"

#include <fstream>
#include <iostream>

#include <boost/program_options.hpp>
#include <boost/type_erasure/any_cast.hpp>

namespace po = boost::program_options;

void CommandLineArguments::parse(int argc, char** argv) {
  po::options_description optionsDescription("Allowed options");
  optionsDescription.add_options()("help,h", "Print usage message.")(
      "FILE,f", po::value<std::string>()->default_value(""),
      "Provide path for input file.")(
      "NUM,n", po::value<unsigned int>()->default_value(500),
      "Number of groups to iterate in seed finding.")(
      "DEVICE,d", po::value<std::string>()->default_value(""),
      "Provide a substring of the preferred device.")(
      "LIST,l", "List available SYCL platforms and devices.")(
      "GPU,G", po::bool_switch(), "Execute code only on gpu. Default is 0.")(
      "ALL,a", po::bool_switch(), "Analyze all groups. Default is 0.")(
      "MATCH,m", po::bool_switch(), "Count seed matches. Default is 0.")(
      "CSV,c", po::bool_switch(), "Output results in csv format");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, optionsDescription), vm);
  po::notify(vm);

  if (vm.count("help") != 0) {
    std::cout << optionsDescription << "\n";
    exit(0);
  }

  if (vm.count("LIST") != 0) {
    Acts::Sycl::listPlatforms();
    exit(0);
  }

  onlyGpu = vm["GPU"].as<bool>();
  matches = vm["MATCH"].as<bool>();
  groups = vm["NUM"].as<unsigned int>();
  deviceName = vm["DEVICE"].as<std::string>();
  allgroup = vm["ALL"].as<bool>();
  csvFormat = vm["CSV"].as<bool>();
  inpFileName = vm["FILE"].as<std::string>();
  std::ifstream s(inpFileName);
  inpFileExists = s.good();
}

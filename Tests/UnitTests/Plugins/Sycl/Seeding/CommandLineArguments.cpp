// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "CommandLineArguments.h"
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/type_erasure/any_cast.hpp>

namespace po = boost::program_options;

void CommandLineArguments::parse(int argc, char** argv) {
    po::options_description optionsDescription("Allowed options");
    optionsDescription.add_options()
      ("help,h","Print usage message.")
      ("inputfile,f", po::value<std::string>(), "Provide path for input file.")
      ("only_gpu,g",po::bool_switch(),"Execute code only on gpu.")
      ("all_groups,a", po::bool_switch(), "Execute on all groups.")
      ("groups,c",po::value<unsigned int>()->default_value(500),"Add number of groups to execute on.")
      ("matches,m",po::bool_switch(),"Count seed matches.")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, optionsDescription), vm);
    po::notify(vm);

    if (vm.count("help") != 0) {  
      std::cout << optionsDescription << "\n";
      exit(0);
    }

    only_gpu = vm["only_gpu"].as<bool>();
    matches = vm["matches"].as<bool>();
    groups = vm["groups"].as<unsigned int>();
    allgroup = vm["all_groups"].as<bool>();
    filename = vm["inputfile"].as<std::string>();
    std::ifstream s(filename);
    fileExists = s.good();
}
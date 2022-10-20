// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/Options.hpp"

#include <exception>
#include <iostream>
#include <string>

#include <TApplication.h>
#include <boost/program_options.hpp>

#define BOOST_AVAILABLE 1
#if ((BOOST_VERSION / 100) % 1000) <= 71
// Boost <=1.71 and lower do not have progress_display.hpp as a replacement yet
#include <boost/progress.hpp>
using progress_display = boost::progress_display;
#else
// Boost >=1.72 can use this as a replacement
#include <boost/timer/progress_display.hpp>
using progress_display = boost::timer::progress_display;
#endif

#include "materialComposition.C"

using namespace boost::program_options;
using VariableReals = ActsExamples::Options::VariableReals;

int main(int argc, char** argv) {
  std::cout << "*** Material Composition plotting " << std::endl;

  try {
    options_description description("*** Usage:");

    // Add the program options
    auto ao = description.add_options();
    ao("help,h", "Display this help message");
    ao("silent,s", bool_switch(), "Silent mode (without X-window/display).");
    ao("input,i", value<std::string>()->default_value(""),
       "Input ROOT file containing the input TTree.");
    ao("tree,t", value<std::string>()->default_value("material-tracks"),
       "Input TTree name.");
    ao("output,o", value<std::string>()->default_value(""),
       "Output ROOT file with histograms");
    ao("bins,b", value<unsigned int>()->default_value(60),
       "Number of bins in eta/phi");
    ao("eta,e", value<float>()->default_value(4.), "Eta range.");
    ao("sub-names", value<std::vector<std::string>>()->multitoken(),
       "Subdetector names.");
    ao("sub-rmin", value<VariableReals>(), "Minimal radial restrictions.");
    ao("sub-rmax", value<VariableReals>(), "Maximal radial restrictions.");
    ao("sub-zmin", value<VariableReals>(), "Minimal z radial restrictions");
    ao("sub-zmax", value<VariableReals>(), "Maximal z radial restrictions.");

    // Set up the variables map
    variables_map vm;
    store(command_line_parser(argc, argv).options(description).run(), vm);
    notify(vm);

    if (vm.count("help") != 0u) {
      std::cout << description;
    }

    // Parse the parameters
    auto iFile = vm["input"].as<std::string>();
    auto iTree = vm["tree"].as<std::string>();
    auto oFile = vm["output"].as<std::string>();

    // Bins & eta range
    unsigned int bins = vm["bins"].as<unsigned int>();
    float eta = vm["eta"].as<float>();

    // Subdetector configurations
    std::vector<Region> dRegion = {};
    auto snames = vm["sub-names"].as<std::vector<std::string>>();
    auto rmins = vm["sub-rmin"].as<VariableReals>().values;
    auto rmaxs = vm["sub-rmax"].as<VariableReals>().values;
    auto zmins = vm["sub-zmin"].as<VariableReals>().values;
    auto zmaxs = vm["sub-zmax"].as<VariableReals>().values;

    size_t subs = snames.size();

    if (subs != rmins.size() or subs != rmaxs.size() or subs != zmins.size() or
        subs != zmaxs.size()) {
      std::cerr << "Configuration problem." << std::endl;
      return 1;
    }

    // Create the regions
    for (unsigned int is = 0; is < subs; ++is) {
      dRegion.push_back(
          {snames[is], rmins[is], rmaxs[is], zmins[is], zmaxs[is]});
    }

    TApplication* tApp =
        vm["silent"].as<bool>()
            ? nullptr
            : new TApplication("ResidualAndPulls", nullptr, nullptr);

    materialComposition(iFile, iTree, oFile, bins, eta, dRegion);

    if (tApp != nullptr) {
      tApp->Run();
    }

  } catch (std::exception& e) {
    std::cerr << e.what() << "\n";
  }

  std::cout << "*** Done." << std::endl;
  return 0;
}

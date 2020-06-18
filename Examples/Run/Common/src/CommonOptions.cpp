// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Options/CommonOptions.hpp"
#include <exception>
#include <fstream>
#include <regex>
#include <system_error>
#include "ACTFW/Utilities/Options.hpp"

using namespace boost::program_options;

boost::program_options::options_description FW::Options::makeDefaultOptions(
    std::string caption) {
  options_description opt(caption);

  opt.add_options()("help,h", "Produce help message");
  opt.add_options()(
      "loglevel,l", value<size_t>()->default_value(2),
      "The output log level. Please set the wished number (0 = VERBOSE, 1 = "
      "DEBUG, 2 = INFO, 3 = WARNING, 4 = ERROR, 5 = FATAL).");
  opt.add_options()(
      "response-file", value<std::string>()->default_value(""),
      "Configuration file (response file) replacing command line options.");

  return opt;
}

void FW::Options::addSequencerOptions(
    boost::program_options::options_description& opt) {
  // sequencer options
  opt.add_options()("events,n", value<size_t>(),
                    "The number of events to process. If not given, all "
                    "available events will be processed.")(
      "skip", value<size_t>()->default_value(0),
      "The number of events to skip")(
      "jobs,j", value<int>()->default_value(-1),
      "Number of parallel jobs, negative for automatic.");
}

void FW::Options::addRandomNumbersOptions(
    boost::program_options::options_description& opt) {
  opt.add_options()("rnd-seed", value<uint64_t>()->default_value(1234567890u),
                    "Random numbers seed.");
}

void FW::Options::addGeometryOptions(
    boost::program_options::options_description& opt) {
  opt.add_options()("geo-surface-loglevel", value<size_t>()->default_value(3),
                    "The outoput log level for the surface building.")(
      "geo-layer-loglevel", value<size_t>()->default_value(3),
      "The output log level for the layer building.")(
      "geo-volume-loglevel", value<size_t>()->default_value(3),
      "The output log level for the volume building.")(
      "geo-detector-volume", value<read_strings>()->default_value({{}}),
      "Sub detectors for the output writing");
}

void FW::Options::addMaterialOptions(
    boost::program_options::options_description& opt) {
  opt.add_options()(
      "mat-input-type", value<std::string>()->default_value("build"),
      "The way material is loaded: 'none', 'build', 'proto', 'file'.")(
      "mat-input-file", value<std::string>()->default_value(""),
      "Name of the material map input file, supported: '.json' or '.root'.")(
      "mat-output-file", value<std::string>()->default_value(""),
      "Name of the material map output file (without extension).")(
      "mat-output-sensitives", value<bool>()->default_value(true),
      "Write material information of sensitive surfaces.")(
      "mat-output-approaches", value<bool>()->default_value(true),
      "Write material information of approach surfaces.")(
      "mat-output-representing", value<bool>()->default_value(true),
      "Write material information of representing surfaces.")(
      "mat-output-boundaries", value<bool>()->default_value(true),
      "Write material information of boundary surfaces.")(
      "mat-output-volumes", value<bool>()->default_value(true),
      "Write material information of volumes.")(
      "mat-output-dense-volumes", value<bool>()->default_value(false),
      "Write material information of dense volumes.")(
      "mat-output-data", value<bool>()->default_value(true),
      "Output the data field(s).")(
      "mat-output-allmaterial", value<bool>()->default_value(false),
      "Add protoMaterial to all surfaces and volume for the mapping.");
}

void FW::Options::addOutputOptions(
    boost::program_options::options_description& opt) {
  // Add specific options for this example
  opt.add_options()("output-dir", value<std::string>()->default_value(""),
                    "Output directory location.")(
      "output-root", value<bool>()->default_value(false),
      "Switch on to write '.root' output file(s).")(
      "output-csv", value<bool>()->default_value(false),
      "Switch on to write '.csv' output file(s).")(
      "output-obj", value<bool>()->default_value(false),
      "Switch on to write '.obj' ouput file(s).")(
      "output-json", value<bool>()->default_value(false),
      "Switch on to write '.json' ouput file(s).")(
      "output-txt", value<bool>()->default_value(false),
      "Switch on to write '.txt' ouput file(s).");
}

void FW::Options::addInputOptions(
    boost::program_options::options_description& opt) {
  // Add specific options for this example
  opt.add_options()("input-dir", value<std::string>()->default_value(""),
                    "Input directory location.")(
      "input-files", value<read_strings>()->multitoken()->default_value({}),
      "Input files, space separated.")("input-root",
                                       value<bool>()->default_value(false),
                                       "Switch on to read '.root' file(s).")(
      "input-csv", value<bool>()->default_value(false),
      "Switch on to read '.csv' file(s).")("input-obj",
                                           value<bool>()->default_value(false),
                                           "Switch on to read '.obj' file(s).")(
      "input-json", value<bool>()->default_value(false),
      "Switch on to read '.json' file(s).");
}

boost::program_options::variables_map FW::Options::parse(
    const boost::program_options::options_description& opt, int argc,
    char* argv[]) noexcept(false) {
  variables_map vm;
  store(command_line_parser(argc, argv).options(opt).run(), vm);
  notify(vm);

  if (vm.count("response-file") and
      not vm["response-file"].template as<std::string>().empty()) {
    // Load the file and tokenize it
    std::ifstream ifs(vm["response-file"].as<std::string>().c_str());
    if (!ifs) {
      throw(std::system_error(std::error_code(),
                              "Could not open response file."));
    }
    // Read the whole file into a string
    std::stringstream ss;
    ss << ifs.rdbuf();
    std::string rString = ss.str();
    std::vector<std::string> args;
    const std::regex rgx("[ \t\r\n\f]");
    std::sregex_token_iterator iter(rString.begin(), rString.end(), rgx, -1);
    std::sregex_token_iterator end;
    for (; iter != end; ++iter) {
      if (std::string(*iter).empty()) {
        continue;
      }
      args.push_back(*iter);
    }
    // Parse the file and store the options
    store(command_line_parser(args).options(opt).run(), vm);
  }

  // Automatically handle help
  if (vm.count("help")) {
    std::cout << opt << std::endl;
    vm.clear();
  }
  return vm;
}

Acts::Logging::Level FW::Options::readLogLevel(
    const boost::program_options::variables_map& vm) {
  return Acts::Logging::Level(vm["loglevel"].as<size_t>());
}

FW::Sequencer::Config FW::Options::readSequencerConfig(
    const boost::program_options::variables_map& vm) {
  Sequencer::Config cfg;
  cfg.skip = vm["skip"].as<size_t>();
  if (not vm["events"].empty()) {
    cfg.events = vm["events"].as<size_t>();
  }
  cfg.logLevel = readLogLevel(vm);
  cfg.numThreads = vm["jobs"].as<int>();
  if (not vm["output-dir"].empty()) {
    cfg.outputDir = vm["output-dir"].as<std::string>();
  }
  return cfg;
}

// Read the random numbers config.
FW::RandomNumbers::Config FW::Options::readRandomNumbersConfig(
    const boost::program_options::variables_map& vm) {
  FW::RandomNumbers::Config cfg;
  cfg.seed = vm["rnd-seed"].as<uint64_t>();
  return cfg;
}

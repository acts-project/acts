// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Options/HepMC3Options.hpp"

#include "ActsExamples/Utilities/Options.hpp"

#include <string>

#include <boost/program_options.hpp>

void ActsExamples::Options::addHepMC3WriterOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();

  opt("hepmc3-stem", value<std::string>()->default_value("events"),
      "The stem of file names of the output");
}

ActsExamples::HepMC3AsciiWriter::Config
ActsExamples::Options::readHepMC3WriterOptions(
    const ActsExamples::Options::Variables& variables) {
  ActsExamples::HepMC3AsciiWriter::Config writerConfig;
  writerConfig.outputDir = variables["output-dir"].as<std::string>();
  writerConfig.outputStem = variables["hepmc3-stem"].as<std::string>();

  return writerConfig;
}

void ActsExamples::Options::addHepMC3ReaderOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();

  opt("hepmc3-stem", value<std::string>()->default_value("events"),
      "The stem of file names of the input");
}

ActsExamples::HepMC3AsciiReader::Config
ActsExamples::Options::readHepMC3ReaderOptions(
    const ActsExamples::Options::Variables& variables) {
  ActsExamples::HepMC3AsciiReader::Config readerConfig;
  readerConfig.inputDir = variables["input-dir"].as<std::string>();
  readerConfig.inputStem = variables["hepmc3-stem"].as<std::string>();

  return readerConfig;
}

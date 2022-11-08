// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Options/JsonOptionsWriter.hpp"

#include <boost/program_options.hpp>

void ActsExamples::Options::addJsonWriterOptions(
    ActsExamples::Options::Description& desc) {
  using namespace boost::program_options;

  desc.add_options()(
      "json-output-precision",
      value<size_t>()->default_value(std::numeric_limits<float>::max_digits10),
      "Floating number output precision.")("json-write-sf-boundaries",
                                           bool_switch(),
                                           "Write tracking boundary surfaces.")(
      "json-write-sf-approach", bool_switch(),
      "Write tracking geometry approach surfaces.")(
      "json-write-sf-layer", bool_switch(),
      "Write tracking geometry layer surfaces.")(
      "json-write-sf-sensitive", bool_switch(),
      "Write tracking geometry sensitive surfaces.")(
      "json-write-sf-names-only", bool_switch(),
      "Write only names of tracking geometry surfaces.");
}

ActsExamples::JsonSurfacesWriter::Config
ActsExamples::Options::readJsonSurfacesWriterConfig(
    const ActsExamples::Options::Variables& vm) {
  ActsExamples::JsonSurfacesWriter::Config cfg;
  if (not vm["output-dir"].empty()) {
    cfg.outputDir = vm["output-dir"].as<std::string>();
  }
  cfg.outputPrecision = vm["json-output-precision"].as<size_t>();
  cfg.writeLayer = vm["json-write-sf-layer"].as<bool>();
  cfg.writeApproach = vm["json-write-sf-approach"].as<bool>();
  cfg.writeSensitive = vm["json-write-sf-sensitive"].as<bool>();
  cfg.writeBoundary = vm["json-write-sf-boundaries"].as<bool>();
  cfg.writeOnlyNames = vm["json-write-sf-names-only"].as<bool>();

  return cfg;
}

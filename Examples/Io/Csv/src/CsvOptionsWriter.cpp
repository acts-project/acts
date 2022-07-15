// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvOptionsWriter.hpp"

#include <limits>

#include <boost/program_options.hpp>
#include <dfe/dfe_io_dsv.hpp>

void ActsExamples::Options::addCsvWriterOptions(
    ActsExamples::Options::Description& desc) {
  using namespace boost::program_options;

  desc.add_options()(
      "csv-output-precision",
      value<size_t>()->default_value(std::numeric_limits<float>::max_digits10),
      "Floating number output precision.")(
      "csv-tg-perevent", bool_switch(), "Write tracking geometry per event.");
}

ActsExamples::CsvParticleWriter::Config
ActsExamples::Options::readCsvParticleWriterConfig(
    const ActsExamples::Options::Variables& vm) {
  ActsExamples::CsvParticleWriter::Config cfg;
  if (not vm["output-dir"].empty()) {
    cfg.outputDir = vm["output-dir"].as<std::string>();
  }
  cfg.outputPrecision = vm["csv-output-precision"].as<size_t>();
  return cfg;
}

ActsExamples::CsvSimHitWriter::Config
ActsExamples::Options::readCsvSimHitWriterConfig(const Variables& vm) {
  ActsExamples::CsvSimHitWriter::Config cfg;
  if (not vm["output-dir"].empty()) {
    cfg.outputDir = vm["output-dir"].as<std::string>();
  }
  return cfg;
}

ActsExamples::CsvPlanarClusterWriter::Config
ActsExamples::Options::readCsvPlanarClusterWriterConfig(
    const ActsExamples::Options::Variables& vm) {
  ActsExamples::CsvPlanarClusterWriter::Config cfg;
  if (not vm["output-dir"].empty()) {
    cfg.outputDir = vm["output-dir"].as<std::string>();
  }
  cfg.outputPrecision = vm["csv-output-precision"].as<size_t>();
  return cfg;
}

ActsExamples::CsvMeasurementWriter::Config
ActsExamples::Options::readCsvMeasurementWriterConfig(
    const ActsExamples::Options::Variables& vm) {
  ActsExamples::CsvMeasurementWriter::Config cfg;
  if (not vm["output-dir"].empty()) {
    cfg.outputDir = vm["output-dir"].as<std::string>();
  }
  // cfg.outputPrecision = vm["csv-output-precision"].as<size_t>();
  return cfg;
}

ActsExamples::CsvTrackingGeometryWriter::Config
ActsExamples::Options::readCsvTrackingGeometryWriterConfig(
    const ActsExamples::Options::Variables& vm) {
  ActsExamples::CsvTrackingGeometryWriter::Config cfg;
  if (not vm["output-dir"].empty()) {
    cfg.outputDir = vm["output-dir"].as<std::string>();
  }
  cfg.outputPrecision = vm["csv-output-precision"].as<size_t>();
  cfg.writePerEvent = (vm.count("csv-tg-perevent") != 0u);
  return cfg;
}

// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"

#include <boost/program_options.hpp>

ActsExamples::CsvParticleReader::Config
ActsExamples::Options::readCsvParticleReaderConfig(const Variables& vm) {
  ActsExamples::CsvParticleReader::Config cfg;
  if (not vm["input-dir"].empty()) {
    cfg.inputDir = vm["input-dir"].as<std::string>();
  }
  return cfg;
}

ActsExamples::CsvSimHitReader::Config
ActsExamples::Options::readCsvSimHitReaderConfig(const Variables& vm) {
  ActsExamples::CsvSimHitReader::Config cfg;
  if (not vm["input-dir"].empty()) {
    cfg.inputDir = vm["input-dir"].as<std::string>();
  }
  return cfg;
}

ActsExamples::CsvPlanarClusterReader::Config
ActsExamples::Options::readCsvPlanarClusterReaderConfig(const Variables& vm) {
  ActsExamples::CsvPlanarClusterReader::Config cfg;
  if (not vm["input-dir"].empty()) {
    cfg.inputDir = vm["input-dir"].as<std::string>();
  }
  return cfg;
}

// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "DigitizationInput.hpp"

#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Utilities/Options.hpp"

ActsExamples::CsvSimHitReader::Config setupSimHitReading(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer) {
  using namespace ActsExamples;

  // Read some standard options
  auto logLevel = Options::readLogLevel(vars);

  // Read truth hits from CSV files
  auto simHitReaderCfg = Options::readCsvSimHitReaderConfig(vars);
  simHitReaderCfg.inputStem = "hits";
  simHitReaderCfg.outputSimHits = "hits";
  sequencer.addReader(
      std::make_shared<CsvSimHitReader>(simHitReaderCfg, logLevel));

  return simHitReaderCfg;
}

ActsExamples::CsvParticleReader::Config setupParticleReading(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer) {
  using namespace ActsExamples;

  // Read some standard options
  auto logLevel = Options::readLogLevel(vars);

  // Read particles (initial states) from CSV files
  auto particleReader = Options::readCsvParticleReaderConfig(vars);
  particleReader.inputStem = "particles_initial";
  particleReader.outputParticles = "particles_initial";
  sequencer.addReader(
      std::make_shared<CsvParticleReader>(particleReader, logLevel));

  return particleReader;
}

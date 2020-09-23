// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootDigitizationWriter.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimIdentifier.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <ios>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

ActsExamples::RootDigitizationWriter::RootDigitizationWriter(
    const ActsExamples::RootDigitizationWriter::Config& cfg,
    Acts::Logging::Level lvl)
    : WriterT(cfg.inputMeasurements, "RootDigitizationWriter", lvl),
      m_cfg(cfg),
      m_outputFile(cfg.rootFile) {
  // inputClusters is already checked by base constructor
  if (m_cfg.inputSimulatedHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
    }
  }
  m_outputFile->cd();
  if (not m_cfg.smearers.empty()) {
    ACTS_DEBUG("Smearers are present, preparing trees.");
    for (const auto& smearer : m_cfg.smearers) {
      auto sFunction = smearer.second;
      std::visit([&](auto&& sm) {}, sFunction);
    }
  }
}

ActsExamples::RootDigitizationWriter::~RootDigitizationWriter() {
  /// Close the file if it's yours
  if (m_cfg.rootFile == nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootDigitizationWriter::endRun() {
  // Write the tree
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootDigitizationWriter::writeT(
    const AlgorithmContext& ctx,
    const ActsExamples::GeometryIdMultimap<
        Acts::FittableMeasurement<ActsExamples::DigitizedHit>>& measurements) {
  // retrieve simulated hits
  const auto& simHits =
      ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimulatedHits);

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  return ActsExamples::ProcessCode::SUCCESS;
}

// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootPropagationSummaryWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Propagator/ConstrainedStep.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Helpers.hpp>

#include <ios>
#include <memory>
#include <ostream>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

namespace ActsExamples {

RootPropagationSummaryWriter::RootPropagationSummaryWriter(
    const RootPropagationSummaryWriter::Config& cfg, Acts::Logging::Level level)
    : WriterT(cfg.inputSummaryCollection, "RootPropagationSummaryWriter",
              level),
      m_cfg(cfg),
      m_outputFile(cfg.rootFile) {
  // An input collection name and tree name must be specified
  if (m_cfg.inputSummaryCollection.empty()) {
    throw std::invalid_argument("Missing input collection");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
    }
  }
  m_outputFile->cd();

  m_outputTree = new TTree(m_cfg.treeName.c_str(),
                           "TTree from RootPropagationSummaryWriter");
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  // Set the branches
  m_outputTree->Branch("event_nr", &m_eventNr);
  m_outputTree->Branch("track_nr", &m_trackNr);

  m_outputTree->Branch("d0", &m_d0);
  m_outputTree->Branch("z0", &m_z0);
  m_outputTree->Branch("phi", &m_phi);
  m_outputTree->Branch("theta", &m_theta);
  m_outputTree->Branch("qOverP", &m_qOverP);
  m_outputTree->Branch("t", &m_t);

  m_outputTree->Branch("eta", &m_eta);
  m_outputTree->Branch("pt", &m_pt);
  m_outputTree->Branch("p", &m_p);

  m_outputTree->Branch("nSteps", &m_nSteps);
  m_outputTree->Branch("nStepTrials", &m_nStepTrials);
  m_outputTree->Branch("pathLength", &m_pathLength);
}

RootPropagationSummaryWriter::~RootPropagationSummaryWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode RootPropagationSummaryWriter::finalize() {
  // Write the tree
  m_outputFile->cd();
  m_outputTree->Write();
  /// Close the file if it's yours
  if (m_cfg.rootFile == nullptr) {
    m_outputFile->Close();
  }

  ACTS_VERBOSE("Wrote particles to tree '" << m_cfg.treeName << "' in '"
                                           << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}

ProcessCode RootPropagationSummaryWriter::writeT(
    const AlgorithmContext& context, const PropagationSummaries& summaries) {
  // Exclusive access to the tree while writing
  std::scoped_lock lock(m_writeMutex);

  // Get the event number
  m_eventNr = static_cast<int>(context.eventNumber);

  // Loop over the step vector of each test propagation in this
  for (const auto& [trackNr, summary] : Acts::enumerate(summaries)) {
    // Set the track number
    m_trackNr = static_cast<int>(trackNr);

    // Set the initial trajectory parameters
    const auto& startParameters = summary.startParameters;
    m_d0 = static_cast<float>(startParameters.parameters()[Acts::eBoundLoc0]);
    m_z0 = static_cast<float>(startParameters.parameters()[Acts::eBoundLoc1]);
    m_phi = static_cast<float>(startParameters.parameters()[Acts::eBoundPhi]);
    m_theta =
        static_cast<float>(startParameters.parameters()[Acts::eBoundTheta]);
    m_qOverP =
        static_cast<float>(startParameters.parameters()[Acts::eBoundQOverP]);
    m_t = static_cast<float>(startParameters.parameters()[Acts::eBoundTime]);

    // Derived initial trajectory parameters
    m_eta = static_cast<float>(
        Acts::VectorHelpers::eta(startParameters.direction()));
    m_pt = static_cast<float>(startParameters.transverseMomentum());
    m_p = static_cast<float>(startParameters.absoluteMomentum());

    // Stepper statistics
    m_nSteps = static_cast<int>(summary.steps.size());
    m_nStepTrials = static_cast<int>(summary.nStepTrials);
    m_pathLength = static_cast<int>(summary.pathLength);

    m_outputTree->Fill();
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

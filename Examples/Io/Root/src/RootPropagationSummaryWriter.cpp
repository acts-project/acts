// This file is part of the Acts project.
//
// Copyright (C) 2017-2024 CERN for the benefit of the Acts project
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

  m_outputTree->Branch("step_volume_id", &m_stepVolumeID);
  m_outputTree->Branch("step_boundary_id", &m_stepBoundaryID);
  m_outputTree->Branch("step_layer_id", &m_stepLayerID);
  m_outputTree->Branch("step_approach_id", &m_stepApproachID);
  m_outputTree->Branch("step_sensitive_id", &m_stepSensitiveID);
  m_outputTree->Branch("step_material", &m_stepMaterial);
  m_outputTree->Branch("step_g_x", &m_stepX);
  m_outputTree->Branch("step_g_y", &m_stepY);
  m_outputTree->Branch("step_g_z", &m_stepZ);
  m_outputTree->Branch("step_d_x", &m_stepDx);
  m_outputTree->Branch("step_d_y", &m_stepDy);
  m_outputTree->Branch("step_d_z", &m_stepDz);
  m_outputTree->Branch("step_type", &m_stepType);
  m_outputTree->Branch("step_acc", &m_stepAcc);
  m_outputTree->Branch("step_act", &m_stepAct);
  m_outputTree->Branch("step_abt", &m_stepAbt);
  m_outputTree->Branch("step_usr", &m_stepUsr);
  m_outputTree->Branch("step_trials", &m_stepTrials);
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
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = context.eventNumber;

  // Initialize the last total trials
  // This is used to calculate the number of trials per step
  std::size_t lastTotalTrials = 0;

  // Loop over the step vector of each test propagation in this
  for (const auto& [trackNr, summary] : Acts::enumerate(summaries)) {
    // Clear the vectors for each collection
    m_stepVolumeID.clear();
    m_stepBoundaryID.clear();
    m_stepLayerID.clear();
    m_stepApproachID.clear();
    m_stepSensitiveID.clear();
    m_stepMaterial.clear();
    m_stepX.clear();
    m_stepY.clear();
    m_stepZ.clear();
    m_stepDx.clear();
    m_stepDy.clear();
    m_stepDz.clear();
    m_stepType.clear();
    m_stepAcc.clear();
    m_stepAct.clear();
    m_stepAbt.clear();
    m_stepUsr.clear();
    m_stepTrials.clear();

    // Set the track number
    m_trackNr = trackNr;

    // Set the initial trajectory parameters
    const auto& startParameters = summary.startParameters;
    m_d0 = startParameters.parameters()[Acts::eBoundLoc0];
    m_z0 = startParameters.parameters()[Acts::eBoundLoc1];
    m_phi = startParameters.parameters()[Acts::eBoundPhi];
    m_theta = startParameters.parameters()[Acts::eBoundTheta];
    m_qOverP = startParameters.parameters()[Acts::eBoundQOverP];
    m_t = startParameters.parameters()[Acts::eBoundTime];

    // Derived initial trajectory parameters
    m_eta = Acts::VectorHelpers::eta(startParameters.direction());
    m_pt = startParameters.transverseMomentum();
    m_p = startParameters.absoluteMomentum();

    // Stepper statistics
    m_nSteps = summary.steps.size();
    m_nStepTrials = summary.nStepTrials;
    m_pathLength = summary.pathLength;

    // Loop over single steps
    for (const auto& step : summary.steps) {
      const auto& geoID = step.geoID;
      m_stepVolumeID.push_back(geoID.volume());
      m_stepBoundaryID.push_back(geoID.boundary());
      m_stepLayerID.push_back(geoID.layer());
      m_stepApproachID.push_back(geoID.approach());
      m_stepSensitiveID.push_back(geoID.sensitive());

      int material = 0;
      if (step.surface) {
        if (step.surface->surfaceMaterial() != nullptr) {
          material = 1;
        }
      }
      m_stepMaterial.push_back(material);

      // kinematic information
      m_stepX.push_back(step.position.x());
      m_stepY.push_back(step.position.y());
      m_stepZ.push_back(step.position.z());
      auto direction = step.momentum.normalized();
      m_stepDx.push_back(direction.x());
      m_stepDy.push_back(direction.y());
      m_stepDz.push_back(direction.z());

      {
        double accuracy = step.stepSize.accuracy();
        double actor = step.stepSize.value(Acts::ConstrainedStep::actor);
        double aborter = step.stepSize.value(Acts::ConstrainedStep::aborter);
        double user = step.stepSize.value(Acts::ConstrainedStep::user);
        double actAbs = std::abs(actor);
        double accAbs = std::abs(accuracy);
        double aboAbs = std::abs(aborter);
        double usrAbs = std::abs(user);

        // todo - fold with direction
        if (actAbs < accAbs && actAbs < aboAbs && actAbs < usrAbs) {
          m_stepType.push_back(0);
        } else if (accAbs < aboAbs && accAbs < usrAbs) {
          m_stepType.push_back(1);
        } else if (aboAbs < usrAbs) {
          m_stepType.push_back(2);
        } else {
          m_stepType.push_back(3);
        }

        // Step size information
        m_stepAcc.push_back(Acts::clampValue<float>(accuracy));
        m_stepAct.push_back(Acts::clampValue<float>(actor));
        m_stepAbt.push_back(Acts::clampValue<float>(aborter));
        m_stepUsr.push_back(Acts::clampValue<float>(user));
      }

      // Stepper efficiency
      m_stepTrials.push_back(step.nTotalTrials - lastTotalTrials);
      lastTotalTrials = step.nTotalTrials;
    }

    m_outputTree->Fill();
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples

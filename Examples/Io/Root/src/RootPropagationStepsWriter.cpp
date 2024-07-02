// This file is part of the Acts project.
//
// Copyright (C) 2017-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootPropagationStepsWriter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
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

ActsExamples::RootPropagationStepsWriter::RootPropagationStepsWriter(
    const ActsExamples::RootPropagationStepsWriter::Config& cfg,
    Acts::Logging::Level level)
    : WriterT(cfg.collection, "RootPropagationStepsWriter", level),
      m_cfg(cfg),
      m_outputFile(cfg.rootFile) {
  // An input collection name and tree name must be specified
  if (m_cfg.collection.empty()) {
    throw std::invalid_argument("Missing input collection");
  } else if (m_cfg.treeName.empty()) {
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
                           "TTree from RootPropagationStepsWriter");
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  // Set the branches
  m_outputTree->Branch("event_nr", &m_eventNr);
  m_outputTree->Branch("volume_id", &m_volumeID);
  m_outputTree->Branch("boundary_id", &m_boundaryID);
  m_outputTree->Branch("layer_id", &m_layerID);
  m_outputTree->Branch("approach_id", &m_approachID);
  m_outputTree->Branch("sensitive_id", &m_sensitiveID);
  m_outputTree->Branch("material", &m_material);
  m_outputTree->Branch("g_x", &m_x);
  m_outputTree->Branch("g_y", &m_y);
  m_outputTree->Branch("g_z", &m_z);
  m_outputTree->Branch("d_x", &m_dx);
  m_outputTree->Branch("d_y", &m_dy);
  m_outputTree->Branch("d_z", &m_dz);
  m_outputTree->Branch("type", &m_step_type);
  m_outputTree->Branch("step_acc", &m_step_acc);
  m_outputTree->Branch("step_act", &m_step_act);
  m_outputTree->Branch("step_abt", &m_step_abt);
  m_outputTree->Branch("step_usr", &m_step_usr);
  m_outputTree->Branch("nStepTrials", &m_nStepTrials);
}

ActsExamples::RootPropagationStepsWriter::~RootPropagationStepsWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootPropagationStepsWriter::finalize() {
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

ActsExamples::ProcessCode ActsExamples::RootPropagationStepsWriter::writeT(
    const AlgorithmContext& context,
    const std::vector<PropagationSteps>& stepCollection) {
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = context.eventNumber;

  // Initialize the last total trials
  // This is used to calculate the number of trials per step
  std::size_t lastTotalTrials = 0;

  // Loop over the step vector of each test propagation in this
  for (auto& steps : stepCollection) {
    // Clear the vectors for each collection
    m_volumeID.clear();
    m_boundaryID.clear();
    m_layerID.clear();
    m_approachID.clear();
    m_sensitiveID.clear();
    m_material.clear();
    m_x.clear();
    m_y.clear();
    m_z.clear();
    m_dx.clear();
    m_dy.clear();
    m_dz.clear();
    m_step_type.clear();
    m_step_acc.clear();
    m_step_act.clear();
    m_step_abt.clear();
    m_step_usr.clear();
    m_nStepTrials.clear();

    // Loop over single steps
    for (auto& step : steps) {
      const auto& geoID = step.geoID;
      m_sensitiveID.push_back(geoID.sensitive());
      m_approachID.push_back(geoID.approach());
      m_layerID.push_back(geoID.layer());
      m_boundaryID.push_back(geoID.boundary());
      m_volumeID.push_back(geoID.volume());

      int material = 0;
      if (step.surface) {
        if (step.surface->surfaceMaterial() != nullptr) {
          material = 1;
        }
      }
      m_material.push_back(material);

      // kinematic information
      m_x.push_back(step.position.x());
      m_y.push_back(step.position.y());
      m_z.push_back(step.position.z());
      auto direction = step.momentum.normalized();
      m_dx.push_back(direction.x());
      m_dy.push_back(direction.y());
      m_dz.push_back(direction.z());

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
        m_step_type.push_back(0);
      } else if (accAbs < aboAbs && accAbs < usrAbs) {
        m_step_type.push_back(1);
      } else if (aboAbs < usrAbs) {
        m_step_type.push_back(2);
      } else {
        m_step_type.push_back(3);
      }

      // Step size information
      m_step_acc.push_back(Acts::clampValue<float>(accuracy));
      m_step_act.push_back(Acts::clampValue<float>(actor));
      m_step_abt.push_back(Acts::clampValue<float>(aborter));
      m_step_usr.push_back(Acts::clampValue<float>(user));

      // Stepper efficiency
      m_nStepTrials.push_back(step.nTotalTrials - lastTotalTrials);
      lastTotalTrials = step.nTotalTrials;
    }
    m_outputTree->Fill();
  }
  return ActsExamples::ProcessCode::SUCCESS;
}

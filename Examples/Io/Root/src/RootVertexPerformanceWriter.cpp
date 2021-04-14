// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootVertexPerformanceWriter.hpp"

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"


#include <ios>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;

ActsExamples::RootVertexPerformanceWriter::RootVertexPerformanceWriter(
    const ActsExamples::RootVertexPerformanceWriter::Config& cfg,
    Acts::Logging::Level lvl)
    : WriterT(cfg.inputVertices, "RootVertexPerformanceWriter", lvl),
      m_cfg(cfg),
      m_outputFile(cfg.rootFile) {

  if (cfg.outputFilename.empty()) {
    throw std::invalid_argument("Missing output filename");
  }
  if (m_cfg.outputTreename.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    auto path = joinPaths(m_cfg.outputDir, m_cfg.outputFilename);
    m_outputFile = TFile::Open(path.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + path);
    }
  }
  m_outputFile->cd();
  m_outputTree =
      new TTree(m_cfg.outputTreename.c_str(), m_cfg.outputTreename.c_str());
  if (m_outputTree == nullptr)
    throw std::bad_alloc();
  else {
    // I/O parameters
    m_outputTree->Branch("diffx", &m_diffx);
    m_outputTree->Branch("diffy", &m_diffy);
    m_outputTree->Branch("diffz", &m_diffz);
    m_outputTree->Branch("recoVtxx", &m_recoVtxx);
    m_outputTree->Branch("recoVtxy", &m_recoVtxy);
    m_outputTree->Branch("recoVtxz", &m_recoVtxz);
    m_outputTree->Branch("trueVtxx", &m_trueVtxx);
    m_outputTree->Branch("trueVtxy", &m_trueVtxy);
    m_outputTree->Branch("trueVtxz", &m_trueVtxz);
    m_outputTree->Branch("nRecoVtx", &m_nrecoVtx);
    m_outputTree->Branch("nTrueVtx", &m_ntrueVtx);
    m_outputTree->Branch("nMaxAcceptanceVtx", &m_nmaxAcceptanceVtx);
  }
}

ActsExamples::RootVertexPerformanceWriter::
    ~RootVertexPerformanceWriter() {
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode
ActsExamples::RootVertexPerformanceWriter::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_outputTree->Write();
  }
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootVertexPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const std::vector<Acts::Vertex<Acts::BoundTrackParameters>>& vertices) {

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  ACTS_DEBUG("Number of reco vertices in event: " << vertices.size());
  if (m_outputFile == nullptr)
    return ProcessCode::SUCCESS;

  // Read input collections
  const auto& truthParticles =
      ctx.eventStore.get<std::vector<ActsFatras::Particle>>(m_cfg.inputTruthParticles);

  std::vector<Acts::Vector3> truePosVec;
  int oldVtxId = -1;
  std::vector<int> allPriVtxIds;
  for(const auto& p : truthParticles){

    int priVtxId = p.particleId().vertexPrimary();
    int secVtxId = p.particleId().vertexSecondary();

    if(secVtxId != 0){
      // truthparticle from secondary vtx
      continue;
    }

    if(priVtxId != oldVtxId){
      const auto& truePos = p.position();
      truePosVec.push_back(truePos);
      oldVtxId = priVtxId;
      allPriVtxIds.push_back(priVtxId);

      m_trueVtxx.push_back(truePos[0]);
      m_trueVtxy.push_back(truePos[1]);
      m_trueVtxz.push_back(truePos[2]);

    }
  }

  std::vector<int> usedIdxs;
  for(const auto& rv : vertices){
    const auto& recoPos = rv.position();

    m_recoVtxx.push_back(recoPos[0]);
    m_recoVtxy.push_back(recoPos[1]);
    m_recoVtxz.push_back(recoPos[2]);

    int clostestIdx = -1;
    double minDist = std::numeric_limits<double>::max();
    for(unsigned int i = 0; i<truePosVec.size(); i++){
      if(std::find(usedIdxs.begin(), usedIdxs.end(), i) != usedIdxs.end()) {
        // True vertex has already been used
        continue;
      } 

      double dist = (recoPos - truePosVec[i]).norm();
      if(dist < minDist){
        clostestIdx = i;
        minDist = dist;
      }
    }
    if(clostestIdx == -1){
      ACTS_WARNING("No matching true vertex found.");
      continue;
    }
    usedIdxs.push_back(clostestIdx);
    m_diffx.push_back(recoPos[0] - truePosVec[clostestIdx][0]);
    m_diffy.push_back(recoPos[1] - truePosVec[clostestIdx][1]);
    m_diffz.push_back(recoPos[2] - truePosVec[clostestIdx][2]);
  }

  m_nrecoVtx = vertices.size();
  m_nmaxAcceptanceVtx = truePosVec.size();
  m_ntrueVtx = *std::max_element(allPriVtxIds.begin(), allPriVtxIds.end());

  ACTS_DEBUG("Number of truth vertices in event: " << m_ntrueVtx);
  ACTS_DEBUG("Number of max. vertices in acceptance in event: " << m_nmaxAcceptanceVtx);
  ACTS_DEBUG("Number of reco vertices in event: " << m_ntrueVtx);

  // fill the variables
  m_outputTree->Fill();

  m_diffx.clear();
  m_diffy.clear();
  m_diffz.clear(); 

  return ProcessCode::SUCCESS;
}

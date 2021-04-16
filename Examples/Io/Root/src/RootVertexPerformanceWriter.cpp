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
  if (m_cfg.inputAllTruthParticles.empty()) {
    throw std::invalid_argument("Collection with all truth particles missing");
  }
  if (m_cfg.inputSelectedTruthParticles.empty()) {
    throw std::invalid_argument("Collection with selected truth particles missing");
  }
  if (m_cfg.inputAssociatedTruthParticles.empty()) {
    throw std::invalid_argument("Collection with track-associated truth particles missing");
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
    m_outputTree->Branch("nVtxDetectorAcceptance", &m_nVtxDetAcceptance);
    m_outputTree->Branch("nVtxReconstructable", &m_nVtxReconstructable);
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

int ActsExamples::RootVertexPerformanceWriter::getNumberOfReconstructableVertices(const SimParticleContainer& collection) const
{
    // map for finding frequency
    std::map<int,int>fmap;
    
    std::vector<int> reconstructableTruthVertices;
   
    // traverse the array for frequency
    for(const auto& p : collection){
      int secVtxId = p.particleId().vertexSecondary();
      if(secVtxId != 0){
       // truthparticle from secondary vtx
        continue;
      }
      int priVtxId = p.particleId().vertexPrimary();
      fmap[priVtxId]++;
    }
   
    // iterate over the map
    for(auto it:fmap)
    {
        // Require at least 2 tracks
        if(it.second > 1)
        {
            reconstructableTruthVertices.push_back(it.first);
        }
    }
    
    return reconstructableTruthVertices.size();
}

int ActsExamples::RootVertexPerformanceWriter::getNumberOfTruePriVertices(const SimParticleContainer& collection) const
{
  // Vector to store indices of all primary vertices
  std::set<int> allPriVtxIds;
  for(const auto& p : collection){

    int priVtxId = p.particleId().vertexPrimary();
    int secVtxId = p.particleId().vertexSecondary();
    if(secVtxId != 0){
      // truthparticle from secondary vtx
      continue;
    }
    // Insert to set, removing duplicates
    allPriVtxIds.insert(priVtxId);
  }
  // Size of set corresponds to total number of primary vertices
  return allPriVtxIds.size();
  
}

ActsExamples::ProcessCode ActsExamples::RootVertexPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const std::vector<Acts::Vertex<Acts::BoundTrackParameters>>& vertices) {

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  m_nrecoVtx = vertices.size();

  ACTS_DEBUG("Number of reco vertices in event: " << m_nrecoVtx);
  if (m_outputFile == nullptr)
    return ProcessCode::SUCCESS;

  // Read truth particle input collection
  const auto& allTruthParticles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputAllTruthParticles);
  // Get number of true primary vertices
  m_ntrueVtx = getNumberOfTruePriVertices(allTruthParticles);

  ACTS_INFO("Total number of generated truth particles in event : " << allTruthParticles.size());
  ACTS_INFO("Total number of generated truth primary vertices : " << m_ntrueVtx);

  // Read selected truth particle input collection
  const auto& allSelectedParticles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputSelectedTruthParticles);
  // Get number of true primary vertices
  m_nVtxDetAcceptance = getNumberOfTruePriVertices(allSelectedParticles);

  ACTS_INFO("Total number of selected truth particles in event : " << allSelectedParticles.size());
  ACTS_INFO("Total number of detector-accepted truth primary vertices : " << m_nVtxDetAcceptance);

  // Read track-associated truth particle input collection
  const auto& allAssociatedTruthParticles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputAssociatedTruthParticles);
  // Get number of true primary vertices
  m_nVtxReconstructable = getNumberOfReconstructableVertices(allAssociatedTruthParticles);

  ACTS_INFO("Total number of reco track-associated truth particles in event : " << allAssociatedTruthParticles.size());
  ACTS_INFO("Total number of reco track-associated truth primary vertices : " << m_nVtxReconstructable);

  // Collect all true primary vertex positions
  std::vector<Acts::Vector3> truePosVec;
  int oldVtxId = -1;
  // Loop over all truth particles associated to
  // fitted tracks
  for(const auto& p : allTruthParticles){

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

      m_trueVtxx.push_back(truePos[0]);
      m_trueVtxy.push_back(truePos[1]);
      m_trueVtxz.push_back(truePos[2]);

    }
  }

  // Delta-r matching between reco and truth vertices
  // Indices of already used truth vertices
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


  // fill the variables
  m_outputTree->Fill();

  m_diffx.clear();
  m_diffy.clear();
  m_diffz.clear(); 

  return ProcessCode::SUCCESS;
}

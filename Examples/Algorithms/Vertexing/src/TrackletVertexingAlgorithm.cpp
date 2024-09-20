// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/TrackletVertexingAlgorithm.hpp"

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/MultiIndex.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Range.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "TH1D.h"
#include "TROOT.h"
#include "TRandom.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include <Eigen/Dense>

using namespace Eigen;

namespace ActsExamples {
struct AlgorithmContext;
}  // namespace ActsExamples

ActsExamples::TrackletVertexingAlgorithm::TrackletVertexingAlgorithm(
    ActsExamples::TrackletVertexingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("TrackletVertexingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {

  if (m_cfg.inputSpacePointsMC.empty()) {
    throw std::invalid_argument("Missing seeds or space point collection");
  }

  for (const auto& spName : m_cfg.inputSpacePointsMC) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputSpacePointsMC.emplace_back(
        std::make_unique<ReadDataHandle<SimSpacePointContainer>>(
            this,
            "InputSpacePoints#" + std::to_string(m_inputSpacePointsMC.size())));
    handle->initialize(spName);
  }

  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_outputRecPrimaryVertex.initialize(m_cfg.outputRecPrimaryVertex);
  m_outputFitPrimaryVertex.initialize(m_cfg.outputFitPrimaryVertex);
  m_outputGenPrimaryVertex.initialize(m_cfg.outputGenPrimaryVertex);
  m_outputFitFunction.initialize(m_cfg.outputFitFunction);
  m_outputZTracklets.initialize(m_cfg.outputZTracklets);
  m_outputZTrackletsPeak.initialize(m_cfg.outputZTrackletsPeak);

  hist = new TH1D("h",";;",m_cfg.nbins,-65,10);

  histMC = new TH1D("hMC",";;",m_cfg.nbins,-65,10);
}

ActsExamples::ProcessCode ActsExamples::TrackletVertexingAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // REAL PART
  const auto& spacePoints = m_inputSpacePoints(ctx);
  const auto& particles = m_inputParticles(ctx);
  hist->Reset("ICESM");
  histMC->Reset("ICESM");

  if(m_cfg.noGuessing){
    double zper = m_cfg.zPerigee;
    m_outputRecPrimaryVertex(ctx, std::move(zper));
    for (const auto& particle : particles) {
      m_outputGenPrimaryVertex(ctx, std::move(particle.position()[2]));
      break;
    }
    return ActsExamples::ProcessCode::SUCCESS;
  }


  double ibz = 71.175; //spacePoints[ib].z()
  double imz = 151.175; //imz
  
  int ev = int(gRandom->Rndm()*10);

  std::vector<std::vector<float>> evhits1;
  std::vector<std::vector<float>> evhits2;
  
  for (size_t ib = 0; ib < spacePoints.size() - 2; ++ib) {
    std::vector<float> sp = {static_cast<float>(spacePoints[ib].x()),static_cast<float>(spacePoints[ib].y())};
    auto z = spacePoints[ib].z();
    if(int(z) == int(ibz)){
      evhits1.push_back(sp);
    }
    
    if(int(z) == int(imz)){
      evhits2.push_back(sp);
    }
  }


  for(auto hits1 : evhits1){
  //for (size_t ib = 0; ib < spacePoints.size() - 2; ++ib) {
    
    if (ibz> m_cfg.zmax)
      continue;//break
    if (ibz < m_cfg.zmin)
      continue;//break
    
    //double r = std::sqrt(spacePoints[ib].x()*spacePoints[ib].x() + spacePoints[ib].y()*spacePoints[ib].y() + ibz*ibz);
    //double phi0 = std::atan2(spacePoints[ib].y(), spacePoints[ib].x());  // atan2 is used to handle the correct quadrant
    double r = std::sqrt(hits1[0]*hits1[0] + hits1[1]*hits1[1] + ibz*ibz);
    double phi0 = std::atan2(hits1[1], hits1[0]);  // atan2 is used to handle the correct quadrant
    double theta0 = std::acos(ibz / r);
    
    //for (size_t im = ib + 1; im < spacePoints.size() - 1; ++im) {
    for(auto hits2 : evhits2){
      if (ibz>=imz)
        continue;
      
      if (imz> m_cfg.zmax)
        continue;//break

      //double phi1 = std::atan2(spacePoints[im].y(), spacePoints[im].x());  // atan2 is used to handle the correct quadrant
      double phi1 = std::atan2(hits2[1], hits2[0]);  // atan2 is used to handle the correct quadrant
       if(abs(phi0-phi1) > m_cfg.deltaPhi)
        continue;

      //r = std::sqrt(spacePoints[im].x()*spacePoints[im].x() + spacePoints[im].y()*spacePoints[im].y() + imz*imz);
      r = std::sqrt(hits2[0]*hits2[0] + hits2[1]*hits2[1] + imz*imz);
      double theta1 = std::acos(imz / r);

      if(theta1-theta0 > m_cfg.deltaThetaMax || theta1-theta0 < m_cfg.deltaThetaMin)
        continue;
        
      double vy = (hits2[1]-hits1[1])/
                (imz-ibz);
      double py = hits1[1]-ibz*vy;

        
      //double vy = (spacePoints[im].y()-spacePoints[ib].y())/
      //          (imz-ibz);
      
      //if(vy*spacePoints[ib].y()<0)
      //  continue;
      //double py = spacePoints[ib].y()-ibz*vy;

      if(m_cfg.projective){
        if(abs(hits2[0]) < 70 && abs(hits2[1]) < 70 && abs(hits1[1]) < 40 && abs(hits1[0]) < 40)
          hist->Fill(-py/vy);
      }
      else{
        hist->Fill(-py/vy);
      }
    }
  }

  // Find the bin with the maximum number of entries
  int maxBin = hist->GetMaximumBin();
  int bin_min =  hist->FindBin(-65);
  int bin_max = hist->FindBin(1);

  // Initialize variables to track the maximum bin
  maxBin = -1;
  double max_value = -1;
  // Loop over the bins in the specified range
  for (int bin = bin_min; bin <= bin_max; ++bin) {
    double value = hist->GetBinContent(bin);
    if (value > max_value) {
        max_value = value;
        maxBin = bin;
    }
  }

  double zPV = hist->GetBinCenter(maxBin);

  std::vector<double> fitParams;
  for(int i=0; i<5; i++)
    fitParams.push_back(0);//fitFunc->GetParameter(i));
  m_outputFitFunction(ctx, std::move(fitParams));

  std::vector<double> zTrackletBins;
  for(int i=0; i<60; i++)
    zTrackletBins.push_back(hist->GetBinContent(i+1));
  m_outputZTracklets(ctx, std::move(zTrackletBins));

  if(m_cfg.doMCtruth){

    //MC part
    // prepare input collections
    const auto& hitParticlesMap = m_inputMeasurementParticlesMap(ctx);
    // compute particle_id -> {hit_id...} map from the
    // hit_id -> {particle_id...} map on the fly.
    const auto& particleHitsMap = invertIndexMultimap(hitParticlesMap);

    // construct the combined input container of space point pointers from all
    // configured input sources.
    // pre-compute the total size required so we only need to allocate once
    size_t nSpacePoints = 0;
    for (const auto& isp : m_inputSpacePointsMC) {
      nSpacePoints += (*isp)(ctx).size();
    }

    std::vector<const SimSpacePoint*> spacePointPtrs;
    spacePointPtrs.reserve(nSpacePoints);
    for (const auto& isp : m_inputSpacePointsMC) {
      for (const auto& spacePoint : (*isp)(ctx)) {
        // since the event store owns the space points, their pointers should be
        // stable and we do not need to create local copies.
        spacePointPtrs.push_back(&spacePoint);
      }
    }

    std::unordered_map<Index, const SimSpacePoint*> spMap;

    for (const auto& spp : spacePointPtrs) {
      if (spp->sourceLinks().empty()) {
        ACTS_WARNING("Missing source link in space point");
        continue;
      }
      for (const auto& slink : spp->sourceLinks()) {
        const IndexSourceLink& islink = slink.get<IndexSourceLink>();
        spMap.emplace(islink.index(), spp);
      }
    }

    for (const auto& particle : particles) {
      m_outputGenPrimaryVertex(ctx, std::move(particle.position()[2]));
      break;
    }

    for (const auto& particle : particles) {
      // find the corresponding hits for this particle
      const auto& hits =
          makeRange(particleHitsMap.equal_range(particle.particleId()));
      // fill hit indices to create the proto track
      ProtoTrack track;
      track.reserve(hits.size());
      for (const auto& hit : hits) {
        track.push_back(hit.second);
      }

      // Space points on the proto track
      std::vector<const SimSpacePoint*> spacePointsOnTrack;
      spacePointsOnTrack.reserve(track.size());
      // Loop over the hit index on the proto track to find the space points
      for (const auto& hitIndex : track) {
        auto it = spMap.find(hitIndex);
        if (it != spMap.end()) {
          spacePointsOnTrack.push_back(it->second);
        }
      }

      if(spacePointsOnTrack.size()<2)
        continue;
      
      std::sort(spacePointsOnTrack.begin(), spacePointsOnTrack.end(),
                [](const SimSpacePoint* lhs, const SimSpacePoint* rhs) {
                  return lhs->z() < rhs->z();
                });
            
      if (spacePointsOnTrack[0]->z()> m_cfg.zmax)
        continue;//break
      if (spacePointsOnTrack[0]->z() < m_cfg.zmin)
        continue;//break  
      if (spacePointsOnTrack[0]->z()>=spacePointsOnTrack[1]->z())
        continue;
      if (spacePointsOnTrack[1]->z()> m_cfg.zmax)
        continue;//break
          
      double vy = (spacePointsOnTrack[1]->y()-spacePointsOnTrack[0]->y())/
                (spacePointsOnTrack[1]->z()-spacePointsOnTrack[0]->z());
      double py = spacePointsOnTrack[0]->y()-spacePointsOnTrack[0]->z()*vy;

      double r = std::sqrt(spacePointsOnTrack[0]->x()*spacePointsOnTrack[0]->x() + spacePointsOnTrack[0]->y()*spacePointsOnTrack[0]->y() + spacePointsOnTrack[0]->z()*spacePointsOnTrack[0]->z());
      double phi0 = std::atan2(spacePointsOnTrack[0]->y(), spacePointsOnTrack[0]->x());  // atan2 is used to handle the correct quadrant
      double theta0 = std::acos(spacePointsOnTrack[0]->z() / r);

      r = std::sqrt(spacePointsOnTrack[1]->x()*spacePointsOnTrack[1]->x() + spacePointsOnTrack[1]->y()*spacePointsOnTrack[1]->y() + spacePointsOnTrack[1]->z()*spacePointsOnTrack[1]->z());
      double phi1 = std::atan2(spacePointsOnTrack[1]->y(), spacePointsOnTrack[1]->x());  // atan2 is used to handle the correct quadrant
      double theta1 = std::acos(spacePointsOnTrack[1]->z() / r);
      if(abs(phi0-phi1) > m_cfg.deltaPhi)
        continue;

      if(theta1-theta0 > m_cfg.deltaThetaMax || theta1-theta0 < m_cfg.deltaThetaMin)
        continue;
      histMC->Fill(-py/vy);
      
    }
    std::vector<double> zTrackletBinsPeak;
    for(int i=0; i<60; i++)
      zTrackletBinsPeak.push_back(histMC->GetBinContent(i+1));
    m_outputZTrackletsPeak(ctx, std::move(zTrackletBinsPeak));
    
  }
  
  m_outputFitPrimaryVertex(ctx, std::move(zPV));
  
  auto zPVId = zPV; //shift in the distribution
  if(zPVId >-7.5)
    zPVId = -0.75;
  else if(zPVId > -21)
    zPVId = -14.25;
  else if(zPVId > -34.5)
    zPVId = -27.75;
  else if(zPVId > -48)
    zPVId = -41.25;
  else
    zPVId = -54.75;
  
  m_outputRecPrimaryVertex(ctx, std::move(zPVId));

  return ActsExamples::ProcessCode::SUCCESS;
}


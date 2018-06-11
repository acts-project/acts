// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include "Acts/Seeding/ICovarianceTool.hpp"
#include "Acts/Seeding/IBinFinder.hpp"
#include "Acts/Seeding/ISeedFilter.hpp"

namespace Acts{

New_Seedmaker::New_Seedmaker(const Acts::Config& config): m_config(config){

}

std::shared_ptr<Acts::SeedmakerState>
New_Seedmaker::initState() const
{
  auto state = std::make_shared<Acts::SeedmakerState>();
  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  state->highland =  13.6*sqrt(m_config.radLengthPerSeed)*(1+0.038*log(m_config.radLengthPerSeed));
  float maxScatteringAngle = state->highland/m_config.minPt;
  state->maxScatteringAngle2 = maxScatteringAngle*maxScatteringAngle;
  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and millimeter
  // TODO: change using ACTS units
  state->pTPerHelixRadius = 300.*m_config.bFieldInZ;
  state->minHelixRadius2 = std::pow(m_config.minPt/state->pTPerHelixRadius,2);
  return state;
}

void 
New_Seedmaker::createSpacePointGrid
( std::vector<const Acts::concept::AnySpacePoint<>*> spVec, 
  std::shared_ptr<Acts::SeedmakerState> state) const
{
  state->seeds.clear();

  SeedingGridConfig gridConf;
  gridConf.bFieldInZ = m_config.bFieldInZ;
  gridConf.minPt = m_config.minPt;
  gridConf.rMax = m_config.rMax;
  gridConf.zMax = m_config.zMax;
  gridConf.zMin = m_config.zMin;
  gridConf.deltaRMax = m_config.deltaRMax;
  gridConf.cotThetaMax = m_config.cotThetaMax;
  // create grid with bin sizes according to the configured geometry
  std::unique_ptr<SPGrid> grid = SPGridCreator::createGrid(gridConf);

  // get region of interest (or full detector if configured accordingly)
  float phiMin = m_config.minPhi;
  float phiMax = m_config.maxPhi;
  float zMin = m_config.zMin;
  float zMax = m_config.zMax;

  // sort by radius
  std::sort(spVec.begin(),spVec.end(),comR());
  for(auto sp : spVec){
    float spZ = sp->z();
    if(spZ > zMax || spZ < zMin) continue;
    float spPhi = std::atan2(sp->y(),sp->x());
    if(spPhi > phiMax || spPhi < phiMin) continue;
    std::array<float,2> cov = m_config.covarianceTool->getCovariances(sp,m_config.zAlign, m_config.rAlign, m_config.sigmaError);
    SPForSeed sps(sp, m_config.beamPos,cov);
    Acts::Vector2D spLocation(spPhi,spZ);
    std::vector<std::shared_ptr<SPForSeed > >& bin = grid->at(spLocation);
    bin.push_back(std::make_shared<SPForSeed>(sps));
  }
  state->binnedSP = std::move(grid);
  state->phiIndex = 1;
  state->zIndex = 1;
}

void
New_Seedmaker::createSeeds
( std::shared_ptr<Acts::SeedmakerState> state) const
{
  std::array<long unsigned int,2ul> phiZbins = state->binnedSP->getNBins();

  
  if(state->outputQueue.size() >= m_config.minSeeds) {
    return;
  }
  bool queueFull = false;

  // loop over all space point bins, break out of loop if queue full
  // store indices in state to continue seed creation after break
  for (; state->phiIndex <= phiZbins[0]; state->phiIndex++){
    for (; state->zIndex <= phiZbins[1]; state->zIndex++){
      std::set<size_t > bottomBins = m_config.bottomBinFinder->findBins(state->phiIndex,state->zIndex,state->binnedSP);
      std::set<size_t > topBins = m_config.topBinFinder->findBins(state->phiIndex,state->zIndex,state->binnedSP);
      createSeedsInRegion(state->binnedSP->at({state->phiIndex,state->zIndex}), bottomBins, topBins, state);
      if(state->outputQueue.size() >= m_config.minSeeds){queueFull = true; break;}
    }
    if(queueFull){break;}
  }
  return;
}

void
New_Seedmaker::createSeedsInRegion
( std::vector<std::shared_ptr<SPForSeed > > currentBin,
  std::set<size_t > bottomBinIndices,
  std::set<size_t > topBinIndices,
  std::shared_ptr<Acts::SeedmakerState> state) const
{

  // loop over all middle space points
  // parallelization requires removal of linCircle* from state
  // and usage of reentrant queue
  // TODO: check how much time reallocation of linCircles in each iteration costs
  for(auto spM : currentBin){
    float rM = spM->radius();
    float zM = spM->z();
    float covrM = spM->covr();
    float covzM = spM->covz();

    std::vector<std::shared_ptr<SPForSeed> > compatBottomSP;

    // bottom space point
    for(auto bottomBinIndex : bottomBinIndices){
      auto bottomBin = state->binnedSP->at(bottomBinIndex);
      for(auto spB : bottomBin){
        float rB = spB->radius();
        float deltaR = rM - rB;
        // if r-distance is too big, try next SP in r-sorted bin
        if (deltaR > m_config.deltaRMax) continue;
        // if r-distance is too small, break because bins are r-sorted
        if (deltaR < m_config.deltaRMin) break;
        // ratio Z/R (forward angle) of space point duplet
        float cotTheta = (zM-spB->z())/deltaR;
        if(std::fabs(cotTheta) > m_config.cotThetaMax) continue;
        // check if duplet origin on z axis within collision region
        float zOrigin = zM-rM*cotTheta;
        if(zOrigin < m_config.collisionRegionMin || zOrigin > m_config.collisionRegionMax) continue;
        compatBottomSP.push_back(spB);
      }
    }
    // no bottom SP found -> try next spM
    if(compatBottomSP.size()==0) continue;
    
    std::vector<std::shared_ptr<SPForSeed> > compatTopSP;

    for(auto topBinIndex : topBinIndices){ 
      auto topBin = state->binnedSP->at(topBinIndex);
      for (auto spT : topBin){
        float rT = spT->radius();
        float deltaR = rT-rM;
        // this condition is the opposite of the condition for bottom SP
        if(deltaR < m_config.deltaRMin ) continue;
        if(deltaR > m_config.deltaRMax ) break;

        float cotTheta = (spT->z()-zM)/deltaR;
        if(std::fabs(cotTheta) > m_config.cotThetaMax) continue;
        float zOrigin = zM-rM*cotTheta;
        if(zOrigin < m_config.collisionRegionMin || zOrigin > m_config.collisionRegionMax) continue;
        compatTopSP.push_back(spT);
      }
    }
    if(compatTopSP.size()==0) continue;
    state->linCircleBottom.clear();
    transformCoordinates(compatBottomSP, spM, true, state->linCircleBottom);
    state->linCircleTop.clear();
    transformCoordinates(compatTopSP, spM, false, state->linCircleTop);
    
    // TODO: significant benefit? avoids compatSp.size()^2 reallocations
    // create vectors here to avoid reallocation in each loop
    std::vector<std::shared_ptr<SPForSeed> > topSpVec;
    std::vector<float > curvatures, impactParameters;

    // TODO: measure cost to reallocate seedsPerSpM each iteration
    std::vector<std::pair<float,std::shared_ptr<InternalSeed > > > seedsPerSpM;

    for(size_t b = 0; b < compatBottomSP.size(); b++){
      auto lb = state->linCircleBottom.at(b);
      float  Zob  = lb.Zo      ;
      float  cotThetaB = lb.cotTheta ;
      float  Vb   = lb.V       ;
      float  Ub   = lb.U       ;
      float  ErB   = lb.Er     ;
      float  iDeltaRB = lb.iDeltaR;

      // 1+(cot^2(theta)) = 1/sin^2(theta)
      float iSinTheta2 = (1.+cotThetaB*cotThetaB) ;
      // calculate max scattering for min momentum at the seed's theta angle
      // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
      // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) < scattering
      // but to avoid trig functions we approximate cot by scaling by 1/sin^4(theta)
      // resolving with pT to p scaling --> only divide by sin^2(theta)
      // max approximation error for allowed scattering angles of 0.04 rad at eta=0: ~8.5%
      float scatteringInRegion2 = state->maxScatteringAngle2 * iSinTheta2;
      // multiply the squared sigma onto the squared scattering
      scatteringInRegion2 *= m_config.sigmaScattering*m_config.sigmaScattering;

      // clear all vectors used in each inner for loop
      topSpVec.clear();
      curvatures.clear();
      impactParameters.clear();
      for(size_t t = 0; t < compatTopSP.size(); t++) {
        auto lt = state->linCircleTop.at(t);

        // add errors of spB-spM and spM-spT pairs and add the correlation term for errors on spM
        float error = lt.Er + ErB + 2*(cotThetaB * lt.cotTheta * covrM + covzM) * iDeltaRB * lt.iDeltaR;

        float deltaCotTheta = cotThetaB - lt.cotTheta;
        float dCotThetaCorrected = deltaCotTheta*deltaCotTheta - error;

        // if deltaTheta larger than the scattering for the lower pT cut, skip
        if ( dCotThetaCorrected > scatteringInRegion2) continue;

        // protects against division by 0
        float dU  = lt.U-Ub; if(dU == 0.) continue ;
        // A and B are evaluated as a function of the circumference parameters x_0 and y_0
        float A   = (lt.V-Vb)/dU                     ;
        float S2  = 1.+A*A                           ;
        float B   = Vb-A*Ub                          ;
        float B2  = B*B                              ;
        // sqrt(S2)/B = 2 * helixradius
        // calculated radius must not be smaller than minimum radius
        if(S2 < B2*state->minHelixRadius2*2) continue;
        // 1/helixradius: (B/sqrt(S2))/2 (we leave everything squared)
        float iHelixradius2 = 4*B2/S2;
        // calculate scattering for p(T) calculated from seed curvature
        float pT2perRadius = state->highland/state->pTPerHelixRadius;
        pT2perRadius = pT2perRadius*pT2perRadius;
        float pT2scatter = iHelixradius2 * pT2perRadius;
        pT2scatter *= pT2scatter;
        //convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta) from rad to deltaCotTheta
        float p2scatter = pT2scatter * iSinTheta2;
        // if deltaTheta larger than allowed scattering for calculated pT, skip
        if(dCotThetaCorrected > p2scatter * m_config.sigmaScattering* m_config.sigmaScattering) continue;
        // A and B allow calculation of impact params in U/V plane with linear function
        // (in contrast to x^2 in x/y plane)
        float Im  = fabs((A-B*rM)*rM)                ;

        if(Im <= m_config.impactMax) {
          // TODO: test impact of at() instead of [] access. should never be out of bounds.
          topSpVec.push_back(compatTopSP.at(t));
          curvatures.push_back(sqrt(iHelixradius2));
          impactParameters.push_back(Im);
        }
      }
      if(!topSpVec.empty()) {
        std::vector<std::pair<float, std::shared_ptr<InternalSeed > > > sameTrackSeeds;
        sameTrackSeeds = m_config.seedFilter->filterSeeds_2SpFixed(compatBottomSP.at(b),
                                                               spM,
                                                               topSpVec,
                                                               curvatures,
                                                               impactParameters,
                                                               Zob);
        seedsPerSpM.insert(seedsPerSpM.end(), sameTrackSeeds.begin(), sameTrackSeeds.end());
      }
    }
    m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, state->outputQueue);
  }
}
  

void New_Seedmaker::transformCoordinates
( std::vector<std::shared_ptr<SPForSeed> >& vec,
  std::shared_ptr<SPForSeed> spM,
  bool bottom,
  std::vector<LinCircle>& linCircleVec) const
{
  float xM = spM->x();
  float yM = spM->y();
  float zM = spM->z();
  float rM = spM->radius();
  float covzM = spM->covz();
  float covrM = spM->covr();
  float cosPhiM = xM/rM;
  float sinPhiM = yM/rM;
  for (auto sp : vec){
    float deltaX = sp->x()-xM;
    float deltaY = sp->y()-yM;
    float deltaZ = sp->z()-zM;
    // calculate projection fraction of spM->sp vector pointing in same direction as
    // vector origin->spM (x) and projection fraction of spM->sp vector pointing 
    // orthogonal to origin->spM (y)
    float x = deltaX * cosPhiM + deltaY*sinPhiM;
    float y = deltaY * cosPhiM - deltaX*sinPhiM;
    // 1/(deltaR*deltaR)
    // x*x+y*y is larger for smaller impact params and higher pT
    // x*x+y*y is always > 1 (unless pT is too low for this Seedfinder)
    // 1/(length of M -> SP) 
    float iDeltaR2 = 1./(deltaX*deltaX+deltaY*deltaY);
    float iDeltaR = sqrt(iDeltaR2);
    // 
    int bottomFactor = 1 * (!bottom) - 1*bottom ;
    // cot_theta = (deltaZ/deltaR)
    float cot_theta = deltaZ*iDeltaR*bottomFactor;
    // VERY frequent (SP^3) access
    LinCircle l;
    l.cotTheta   = cot_theta                                        ;
    // location on z-axis of this SP-duplet
    l.Zo         = zM-rM * cot_theta                                ;
    l.iDeltaR    = iDeltaR                                          ;
    // transformation of circle equation (x,y) into linear equation (u,v)
    // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
    // is transformed into
    // 1 - 2x_0*u - 2y_0*v = 0
    // using the following m_U and m_V
    // (u = A + B*v); A and B are created later on
    l.U    = x*iDeltaR2                                             ;
    l.V    = y*iDeltaR2                                             ;
    //error term for sp-pair without correlation of middle space point
    l.Er   = ((covzM+sp->covz())+(cot_theta*cot_theta)*(covrM+sp->covr()))*iDeltaR2;
    linCircleVec.push_back(l);
  }
}
} // Acts namespace

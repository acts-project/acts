// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include <numeric>

#include "Acts/Seeding/IBinFinder.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include <iostream>

namespace Acts{

New_Seedmaker::New_Seedmaker(const Acts::SeedmakerConfig& config): m_config(config){

  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  m_config.highland =  13.6*sqrt(m_config.radLengthPerSeed)*(1+0.038*log(m_config.radLengthPerSeed));
  float maxScatteringAngle = m_config.highland/m_config.minPt;
  m_config.maxScatteringAngle2 = maxScatteringAngle*maxScatteringAngle;
  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and millimeter
  // TODO: change using ACTS units
  m_config.pTPerHelixRadius = 300.*m_config.bFieldInZ;
  m_config.minHelixRadius2 = std::pow(m_config.minPt/m_config.pTPerHelixRadius,2);
}

template <typename Seed, typename SpacePoint>

std::unique_ptr<Seed> New_Seedmaker::nextSeed(Acts::SeedmakerState* state,
                             std::vector<const SpacePoint*>&     inputSP){
  auto seed = std::make_unique<Seed>();
  if(state->outputQueue.size()>0){
    auto intSeed = std::move(state->outputQueue.front());
    state->outputQueue.pop();
    seed->erase();
    seed->add(inputSP[intSeed->spacepoint0()->spIndex()]);
    seed->add(inputSP[intSeed->spacepoint1()->spIndex()]);
    seed->add(inputSP[intSeed->spacepoint2()->spIndex()]);
    seed->setZVertex(intSeed->z());
  }
  return seed;
}


template <typename SpacePoint>
std::shared_ptr<Acts::SeedmakerState>
New_Seedmaker::initState
(const std::vector<const SpacePoint* >& spVec,
 std::function<Acts::Vector2D(const SpacePoint*,float,float,float)> covTool) const
{
  auto state = std::make_shared<Acts::SeedmakerState>();
  // setup spacepoint grid config
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
  
  std::vector<size_t> indexVec(spVec.size());
  std::iota(std::begin(indexVec),std::end(indexVec),0);

  // sort by radius
  float inverseRBinSize = m_config.deltaRMin * 3;
  // add magnitude of beamPos to rMax to avoid excluding measurements
  size_t numRBins = (m_config.rMax + m_config.beamPos.norm()) * inverseRBinSize;
  std::vector<std::vector<std::unique_ptr<const InternalSpacePoint> > >rBins(numRBins);
  for (size_t spIndex : indexVec){
    const SpacePoint* sp = spVec[spIndex];
    float spX = sp->x();
    float spY = sp->y();
    float spZ = sp->z();

    if(spZ > zMax || spZ < zMin) continue;
    float spPhi = std::atan2(spY,spX);
    if(spPhi > phiMax || spPhi < phiMin) continue;

    // covariance tool provided by user
    Acts::Vector2D cov = covTool(sp,m_config.zAlign, m_config.rAlign, m_config.sigmaError);
    Acts::Vector3D spPosition(spX,spY,spZ);
    auto isp = std::make_unique<const InternalSpacePoint>(spIndex, spPosition, m_config.beamPos, cov);
    // calculate r-Bin index and protect against overflow (underflow not possible)
    int rIndex = isp->radius()*inverseRBinSize;
    // if index out of bounds, the SP is outside the region of interest
    if (rIndex >= numRBins){
      continue;
    }
    rBins[rIndex].push_back(std::move(isp));
  }
//  too expensive :(
//  std::sort(indexVec.begin(),indexVec.end(),[spVec]
//       (const size_t indA, const size_t indB)
//        {auto spA = spVec[indA];
//        auto spB = spVec[indB];
//        // compare x*x+y*y (i.e. r^2) of both SP
//        float rA = spA->x()*spA->x()+spA->y()*spA->y();
//        float rB = spB->x()*spB->x()+spB->y()*spB->y();
//        return rA < rB; });
  for(auto& rbin : rBins){
    for(auto& isp : rbin){
//    done during r-binning
//    for(auto spIndex : indexVec){
//      auto sp = spVec[spIndex];
//      float spX = sp->x();
//      float spY = sp->y();
//      float spZ = sp->z();
//      if(spZ > zMax || spZ < zMin) continue;
//      float spPhi = std::atan2(spY,spX);
//      if(spPhi > phiMax || spPhi < phiMin) continue;
//      // covariance configuration should be done outside of the Seedmaker
//      Acts::Vector2D cov = covTool(sp,m_config.zAlign, m_config.rAlign, m_config.sigmaError);
//      Acts::Vector3D globalPos(spX,spY,spZ);
//      InternalSpacePoint sps(spIndex, globalPos, m_config.beamPos, cov);
      Acts::Vector2D spLocation(isp->phi(),isp->z());
      std::vector<std::unique_ptr<const InternalSpacePoint> >& bin = grid->at(spLocation);
      bin.push_back(std::move(isp));
    }
  }
  state->binnedSP = std::move(grid);
  state->phiIndex = 1;
  state->zIndex = 1;
  return state;
}

void
New_Seedmaker::createSeeds
( std::shared_ptr<Acts::SeedmakerState> state) const
{
  std::array<long unsigned int,2ul> phiZbins = state->binnedSP->getNBins();

  
  if(state->outputQueue.size() >= m_config.minSeeds) {
    return;
  }
  // loop over all space point bins, break out of loop if queue full
  // store indices in state to continue seed creation after break
  for (; state->phiIndex <= phiZbins[0]; state->phiIndex++){
    for (; state->zIndex <= phiZbins[1]; state->zIndex++){
      std::set<size_t > bottomBins = m_config.bottomBinFinder->findBins(state->phiIndex,state->zIndex,state->binnedSP);
      std::set<size_t > topBins = m_config.topBinFinder->findBins(state->phiIndex,state->zIndex,state->binnedSP);
      auto& curbin = state->binnedSP->at({state->phiIndex,state->zIndex});
      createSeedsInRegion(curbin, bottomBins, topBins, state);
      if(state->outputQueue.size() >= m_config.minSeeds){return;}
    }
    state->zIndex = 1;
  }
  return;
}

void
New_Seedmaker::createSeedsInRegion
( std::vector<std::unique_ptr<const InternalSpacePoint > >& currentBin,
  std::set<size_t > bottomBinIndices,
  std::set<size_t > topBinIndices,
  std::shared_ptr<Acts::SeedmakerState> state) const
{

  // loop over all middle space points
  // parallelization requires removal of linCircle* from state
  // and usage of reentrant queue
  // TODO: check how much time reallocation of linCircles in each iteration costs
  for(auto& spMunique : currentBin){
    const InternalSpacePoint* spM = spMunique.get();
    float rM = spM->radius();
    float zM = spM->z();
    float covrM = spM->covr();
    float covzM = spM->covz();

    std::vector<const InternalSpacePoint* > compatBottomSP;

    // bottom space point
    for(auto bottomBinIndex : bottomBinIndices){
      auto& bottomBin = state->binnedSP->at(bottomBinIndex);
      for(auto& spB : bottomBin){
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
        compatBottomSP.push_back(spB.get());
      }
    }
    // no bottom SP found -> try next spM
    if(compatBottomSP.size()==0) continue;
    
    std::vector<const InternalSpacePoint* > compatTopSP;

    for(auto topBinIndex : topBinIndices){ 
      auto& topBin = state->binnedSP->at(topBinIndex);
      for (auto& spT : topBin){
        float rT = spT->radius();
        float deltaR = rT-rM;
        // this condition is the opposite of the condition for bottom SP
        if(deltaR < m_config.deltaRMin ) continue;
        if(deltaR > m_config.deltaRMax ) break;

        float cotTheta = (spT->z()-zM)/deltaR;
        if(std::fabs(cotTheta) > m_config.cotThetaMax) continue;
        float zOrigin = zM-rM*cotTheta;
        if(zOrigin < m_config.collisionRegionMin || zOrigin > m_config.collisionRegionMax) continue;
        compatTopSP.push_back(spT.get());
      }
    }
    if(compatTopSP.size()==0) continue;
    state->linCircleBottom.clear();
    transformCoordinates(compatBottomSP, spM, true, state->linCircleBottom);
    state->linCircleTop.clear();
    transformCoordinates(compatTopSP, spM, false, state->linCircleTop);
    
    // TODO: significant benefit? avoids compatSp.size()^2 reallocations
    // create vectors here to avoid reallocation in each loop
    std::vector<const InternalSpacePoint* > topSpVec;
    std::vector<float > curvatures;
    std::vector<float > impactParameters;

    // TODO: measure cost to reallocate seedsPerSpM each iteration
    std::vector<std::pair<float,std::unique_ptr<const InternalSeed> > > seedsPerSpM;

    for(size_t b = 0; b < compatBottomSP.size(); b++){
      auto lb = state->linCircleBottom[b];
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
      // max approximation error for allowed scattering angles of 0.04 rad at eta=infinity: ~8.5%
      float scatteringInRegion2 = m_config.maxScatteringAngle2 * iSinTheta2;
      // multiply the squared sigma onto the squared scattering
      scatteringInRegion2 *= m_config.sigmaScattering*m_config.sigmaScattering;

      // clear all vectors used in each inner for loop
      topSpVec.clear();
      curvatures.clear();
      impactParameters.clear();
      for(size_t t = 0; t < compatTopSP.size(); t++) {
        auto lt = state->linCircleTop[t];

        // add errors of spB-spM and spM-spT pairs and add the correlation term for errors on spM
        float error2 = lt.Er + ErB + 2*(cotThetaB * lt.cotTheta * covrM + covzM) * iDeltaRB * lt.iDeltaR;

        float deltaCotTheta = std::abs(cotThetaB - lt.cotTheta);

        // if deltaTheta larger than the scattering for the lower pT cut, skip
        float deltaCotTheta2 = deltaCotTheta*deltaCotTheta;
        float error = std::sqrt(error2);
        float dCotThetaMinusError2 = deltaCotTheta2 + error2 - 2*deltaCotTheta*error;
        // avoid taking root of scatteringInRegion
        // if left side of ">" is positive, both sides of unequality can be squared 
        // (scattering is always positive)
        if (deltaCotTheta - error > 0 && dCotThetaMinusError2 > scatteringInRegion2 ) continue;
        //TODO: did everything squared, taking root now. slows down by how much?
        //float dCotThetaCorrected = deltaCotTheta*deltaCotTheta - error;
        //if ( dCotThetaCorrected > scatteringInRegion2) continue;

        // protects against division by 0
        float dU  = lt.U-Ub; if(dU == 0.) continue ;
        // A and B are evaluated as a function of the circumference parameters x_0 and y_0
        float A   = (lt.V-Vb)/dU                     ;
        float S2  = 1.+A*A                           ;
        float B   = Vb-A*Ub                          ;
        float B2  = B*B                              ;
        // sqrt(S2)/B = 2 * helixradius
        // calculated radius must not be smaller than minimum radius
        if(S2 < B2*m_config.minHelixRadius2*4) continue;
        // 1/helixradius: (B/sqrt(S2))/2 (we leave everything squared)
        float iHelixradius2 = 4*B2/S2;
        // calculate scattering for p(T) calculated from seed curvature
        float pT2perRadius = m_config.highland/m_config.pTPerHelixRadius;
        pT2perRadius = pT2perRadius*pT2perRadius;
        float pT2scatter = iHelixradius2 * pT2perRadius;
        //FIXME: include upper pT limit for scatter calc
        //convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta) from rad to deltaCotTheta
        float p2scatter = pT2scatter * iSinTheta2;
        // if deltaTheta larger than allowed scattering for calculated pT, skip
        if(deltaCotTheta - error > 0 && dCotThetaMinusError2 > p2scatter * m_config.sigmaScattering* m_config.sigmaScattering) continue;
        // A and B allow calculation of impact params in U/V plane with linear function
        // (in contrast to x^2 in x/y plane)
        float Im  = std::abs((A-B*rM)*rM)                ;

        if(Im <= m_config.impactMax) {
          topSpVec.push_back(compatTopSP[t]);
          curvatures.push_back(sqrt(iHelixradius2));
          impactParameters.push_back(Im);
        }
      }
      if(!topSpVec.empty()) {
        std::vector<std::pair<float, std::unique_ptr<const InternalSeed> > > sameTrackSeeds;
        sameTrackSeeds = std::move(m_config.seedFilter->filterSeeds_2SpFixed(compatBottomSP[b],
                                                               spM,
                                                               topSpVec,
                                                               curvatures,
                                                               impactParameters,
                                                               Zob));
        seedsPerSpM.insert(seedsPerSpM.end(), std::make_move_iterator(sameTrackSeeds.begin()), std::make_move_iterator(sameTrackSeeds.end()));
      }
    }
    m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, state->outputQueue);
  }
}
  

void New_Seedmaker::transformCoordinates
( std::vector<const InternalSpacePoint* >& vec,
  const InternalSpacePoint * spM,
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

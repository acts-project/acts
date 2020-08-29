// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Sycl/Seeding/Seedfinder.hpp"
#include "Acts/Plugins/Sycl/Seeding/detail/Types.h"
#include <CL/sycl.hpp>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>
#include <array>
#include <exception>
#include <algorithm>

namespace Acts::Sycl {
  /// Kernel classes in order of execution.
  class duplet_search_bottom_kernel;
  class duplet_search_top_kernel;
  class ind_copy_bottom_kernel;
  class ind_copy_top_kernel;
  class transform_coord_bottom_kernel;
  class transform_coord_top_kernel;
  class triplet_search_kernel;
  class filter_2sp_fixed_kernel;

  void offloadComputations(cl::sycl::queue* q,
                          const detail::deviceSeedfinderConfig& configData,
                          const DeviceExperimentCuts& deviceCuts,
                          const std::vector<detail::deviceSpacePoint>& bottomSPs,
                          const std::vector<detail::deviceSpacePoint>& middleSPs,
                          const std::vector<detail::deviceSpacePoint>& topSPs,
                          std::vector<std::vector<SeedData>>& seeds)
  {
    // Each vector stores data of space points in simplified 
    // structures of float variables
    // M: number of middle space points
    // B: number of bottom space points
    // T: number of top space points
    const uint32_t M = middleSPs.size();
    const uint32_t B = bottomSPs.size();
    const uint32_t T = topSPs.size();

    // Up to the Nth space point, the sum of compatible bottom/top space points.
    // We need these for indexing other vectors later in the algorithm.
    std::vector<uint32_t> sumBotCompUptoMid(M+1,0);
    std::vector<uint32_t> sumTopCompUptoMid(M+1,0);
    std::vector<uint32_t> sumBotTopCombined(M+1,0);

    // Similarly to indBotCompMid and indTopCompMid, we store the indices of the 
    // middle space points of the corresponding edges.
    std::vector<uint32_t> indMidBotComp;
    std::vector<uint32_t> indMidTopComp;

    // Number of edges for middle-bottom and middle-top duplet bipartite graphs.
    size_t edgesBottom = 0;
    size_t edgesTop = 0;
    size_t edgesComb = 0;

    try {
      using am = cl::sycl::access::mode;
      using at = cl::sycl::access::target;

      unsigned long globalBufferSize = q->get_device().get_info<cl::sycl::info::device::global_mem_size>();
      unsigned long maxWorkItemPerGroup = q->get_device().get_info<cl::sycl::info::device::max_work_group_size>();  

      // Device allocations
      detail::deviceSpacePoint* botSPArray = cl::sycl::malloc_device<detail::deviceSpacePoint>(B, *q);
      detail::deviceSpacePoint* midSPArray = cl::sycl::malloc_device<detail::deviceSpacePoint>(M, *q);
      detail::deviceSpacePoint* topSPArray = cl::sycl::malloc_device<detail::deviceSpacePoint>(T, *q);

      // Store the number of compatible bottom/top space points per middle space point.
      std::atomic_uint* numBotArray = cl::sycl::malloc_device<std::atomic_uint>(M, *q);
      std::atomic_uint* numTopArray = cl::sycl::malloc_device<std::atomic_uint>(M, *q);

      // The limit of compatible bottom [top] space points per middle space point is B [T].
      // Temporarily we reserve buffers of this size (M*B and M*T).
      // We store the indices of bottom [top] space points in bottomSPs [topSPs].
      // We move the indices to optimal size vectors for easier indexing.
      uint32_t* tmpIndBotArray = cl::sycl::malloc_device<uint32_t>(M*B, *q);
      uint32_t* tmpIndTopArray = cl::sycl::malloc_device<uint32_t>(M*T, *q);

      q->memcpy(botSPArray, bottomSPs.data(), sizeof(detail::deviceSpacePoint)*(B));
      q->memcpy(midSPArray, middleSPs.data(), sizeof(detail::deviceSpacePoint)*(M));
      q->memcpy(topSPArray, topSPs.data(), sizeof(detail::deviceSpacePoint)*(T));
      q->memset(numBotArray, 0, M * sizeof(std::atomic_uint));
      q->memset(numTopArray, 0, M * sizeof(std::atomic_uint));

      q->wait();

      //*********************************************//
      // ********** DUPLET SEARCH - BEGIN ********** //
      //*********************************************//

      // After completing the duplet search, we'll have successfully contructed two
      // bipartite graphs for bottom-middle and top-middle space points.
      // We store the indices of the bottom/top space points of the edges of the graphs.
      // They index the bottomSPs and topSPs vectors.

      {
        q->submit([&](cl::sycl::handler &h) {
          h.parallel_for<duplet_search_bottom_kernel>
            (cl::sycl::range<2>{M,B}, [=](cl::sycl::id<2> idx) {
            const auto mid = idx[0];
            const auto bot = idx[1];
            // We check whether this thread actually makes sense (within bounds).
            if(mid < M && bot < B){

              const auto midSP = midSPArray[mid];
              const auto botSP = botSPArray[bot];

              const auto deltaR = midSP.r - botSP.r;
              const auto cotTheta = (midSP.z - botSP.z) / deltaR;
              const auto zOrigin = midSP.z - midSP.r * cotTheta;

              if( !(deltaR < configData.deltaRMin) &&
                  !(deltaR > configData.deltaRMax) &&
                  !(cl::sycl::abs(cotTheta) > configData.cotThetaMax) &&
                  !(zOrigin < configData.collisionRegionMin) &&
                  !(zOrigin > configData.collisionRegionMax)) {
                // We keep counting duplets with atomic variables 
                const auto ind = numBotArray[mid].fetch_add(1);
                tmpIndBotArray[mid * B + ind] = bot;
              }
            }
          });
        });
    
        q->submit([&](cl::sycl::handler &h) {        
          h.parallel_for<duplet_search_top_kernel>
            (cl::sycl::range<2>{M,B}, [=](cl::sycl::id<2> idx) {
            const auto mid = idx[0];
            const auto top = idx[1];
            // We check whether this thread actually makes sense (within bounds).
            if(mid < M && top < T){

              const auto midSP = midSPArray[mid];
              const auto topSP = topSPArray[top];

              const auto deltaR = topSP.r - midSP.r;
              const auto cotTheta = (topSP.z - midSP.z) / deltaR;
              const auto zOrigin = midSP.z - midSP.r * cotTheta;

              if( !(deltaR < configData.deltaRMin) &&
                  !(deltaR > configData.deltaRMax) &&
                  !(cl::sycl::abs(cotTheta) > configData.cotThetaMax) &&
                  !(zOrigin < configData.collisionRegionMin) &&
                  !(zOrigin > configData.collisionRegionMax)) {
                // We keep counting duplets with atomic access. 
                const auto ind = numTopArray[mid].fetch_add(1);
                tmpIndTopArray[mid * T + ind] = top;
              }
            }
          });
        });
      } // sync

      //*********************************************//
      // *********** DUPLET SEARCH - END *********** //
      //*********************************************//

      // retrieve results from counting duplets
      {
        uint nB[M];
        uint nT[M];

        q->memcpy(&nB[0], numBotArray, M*sizeof(std::atomic_uint));
        q->memcpy(&nT[0], numTopArray, M*sizeof(std::atomic_uint));
        q->wait();

        for(uint32_t i = 1; i < M + 1; ++i){
          sumBotCompUptoMid[i] += sumBotCompUptoMid[i-1] + nB[i-1];
          sumTopCompUptoMid[i] += sumTopCompUptoMid[i-1] + nT[i-1];
          sumBotTopCombined[i] += sumBotTopCombined[i-1] + nT[i-1]*nB[i-1];
        }

        edgesBottom = sumBotCompUptoMid[M];
        edgesTop = sumTopCompUptoMid[M];
        edgesComb = sumBotTopCombined[M];

        if(edgesBottom == 0 || edgesTop == 0) return;

        indMidBotComp.reserve(edgesBottom);
        indMidTopComp.reserve(edgesTop);

        for(uint32_t mid = 0; mid < M; ++mid) {
          std::fill_n(std::back_inserter(indMidBotComp), nB[mid], mid);
          std::fill_n(std::back_inserter(indMidTopComp), nT[mid], mid);
        }
      } // sync

      const auto deviceEdgesBottom = ((edgesBottom/maxWorkItemPerGroup)+1)*maxWorkItemPerGroup;
      const auto deviceEdgesTop =   ((edgesTop/maxWorkItemPerGroup)+1)*maxWorkItemPerGroup;

      const auto edgesBottomOffset = deviceEdgesBottom-edgesBottom;
      const auto edgesTopOffset = deviceEdgesTop-edgesTop;

      cl::sycl::nd_range<1> edgesBotNdRange{
        cl::sycl::range<1>(deviceEdgesBottom),
        cl::sycl::range<1>(maxWorkItemPerGroup)};

      cl::sycl::nd_range<1> edgesTopNdRange{
        cl::sycl::range<1>(deviceEdgesTop),
        cl::sycl::range<1>(maxWorkItemPerGroup)};

      // Copy indices from temporary matrices to final, optimal size vectors.
      uint32_t* indBotArray = cl::sycl::malloc_device<uint32_t>(deviceEdgesBottom,*q);
      uint32_t* indTopArray = cl::sycl::malloc_device<uint32_t>(deviceEdgesTop,*q);
      uint32_t* indMidBotArray = cl::sycl::malloc_device<uint32_t>(deviceEdgesBottom,*q);
      uint32_t* indMidTopArray = cl::sycl::malloc_device<uint32_t>(deviceEdgesTop,*q);
      uint32_t* sumBotArray = cl::sycl::malloc_device<uint32_t>(M+1,*q);
      uint32_t* sumTopArray = cl::sycl::malloc_device<uint32_t>(M+1,*q);
      uint32_t* sumCombArray = cl::sycl::malloc_device<uint32_t>(M+1,*q);
      deviceLinEqCircle* linBotArray = cl::sycl::malloc_device<deviceLinEqCircle>(deviceEdgesBottom,*q);
      deviceLinEqCircle* linTopArray = cl::sycl::malloc_device<deviceLinEqCircle>(deviceEdgesTop, *q); 
      
      q->memset(indBotArray+edgesBottom, 0, edgesBottomOffset * sizeof(uint32_t));
      q->memset(indTopArray+edgesTop, 0, edgesTopOffset * sizeof(uint32_t));
      q->memset(indMidBotArray+edgesBottom, 0, edgesBottomOffset * sizeof(uint32_t));
      q->memset(indMidTopArray+edgesTop, 0, edgesTopOffset * sizeof(uint32_t));

      q->memcpy(indMidBotArray, indMidBotComp.data(), sizeof(uint32_t)*edgesBottom);
      q->memcpy(indMidTopArray, indMidTopComp.data(), sizeof(uint32_t)*edgesTop);
      q->memcpy(sumBotArray, sumBotCompUptoMid.data(), sizeof(uint32_t)*(M+1));
      q->memcpy(sumTopArray, sumTopCompUptoMid.data(), sizeof(uint32_t)*(M+1));
      q->memcpy(sumCombArray, sumBotTopCombined.data(), sizeof(uint32_t)*(M+1));
      
      q->wait();

      {
        q->submit([&](cl::sycl::handler &h){
          h.parallel_for<ind_copy_bottom_kernel>
          (edgesBotNdRange, [=](cl::sycl::nd_item<1> item) {
          auto idx = item.get_global_linear_id();
            uint32_t mid = indMidBotArray[idx];
            uint32_t ind = tmpIndBotArray[mid*B + idx - sumBotArray[mid]];
            indBotArray[idx] = ind;         
          });
        });

        q->submit([&](cl::sycl::handler &h){
          h.parallel_for<ind_copy_top_kernel>
          (edgesTopNdRange, [=](cl::sycl::nd_item<1> item) {
          auto idx = item.get_global_linear_id();
            uint32_t mid = indMidTopArray[idx];
            uint32_t ind = tmpIndTopArray[mid*T + idx - sumTopArray[mid]];
            indTopArray[idx] = ind;     
          });
        });
      } // sync
    
      //************************************************//
      // *** LINEAR EQUATION TRANSFORMATION - BEGIN *** //
      //************************************************//

      // transformation of circle equation (x,y) into linear equation (u,v)
      // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
      // is transformed into
      // 1 - 2x_0*u - 2y_0*v = 0

      // coordinate transformation middle-bottom pairs
      auto linB = q->submit([&](cl::sycl::handler &h) {
        h.parallel_for<transform_coord_bottom_kernel>
          (edgesBotNdRange, [=](cl::sycl::nd_item<1> item) {
          auto idx = item.get_global_linear_id();

          const auto midSP = midSPArray[indMidBotArray[idx]];
          const auto botSP = botSPArray[indBotArray[idx]];

          const auto xM =          midSP.x;
          const auto yM =          midSP.y;
          const auto zM =          midSP.z;
          const auto rM =          midSP.r;
          const auto varianceZM =  midSP.varZ;
          const auto varianceRM =  midSP.varR;
          const auto cosPhiM =     xM / rM;
          const auto sinPhiM =     yM / rM;

          const auto deltaX = botSP.x - xM;
          const auto deltaY = botSP.y - yM;
          const auto deltaZ = botSP.z - zM;

          const auto x = deltaX * cosPhiM + deltaY * sinPhiM;
          const auto y = deltaY * cosPhiM - deltaX * sinPhiM;
          const auto iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
          const auto iDeltaR = cl::sycl::sqrt(iDeltaR2);
          const auto cot_theta = -(deltaZ * iDeltaR);

          deviceLinEqCircle L;
          L.cotTheta = cot_theta;
          L.zo = zM - rM * cot_theta;
          L.iDeltaR = iDeltaR;
          L.u = x * iDeltaR2;
          L.v = y * iDeltaR2;
          L.er = ((varianceZM + botSP.varZ) +
          (cot_theta * cot_theta) * (varianceRM + botSP.varR)) * iDeltaR2;

          linBotArray[idx] = L;
        });
      });
      
      // coordinate transformation middle-top pairs
      auto linT = q->submit([&](cl::sycl::handler &h) {
        h.parallel_for<transform_coord_top_kernel>
          (edgesTopNdRange, [=](cl::sycl::nd_item<1> item) {
          auto idx = item.get_global_linear_id();

          const auto midSP = midSPArray[indMidTopArray[idx]];
          const auto topSP = topSPArray[indTopArray[idx]];

          const auto xM =          midSP.x;
          const auto yM =          midSP.y;
          const auto zM =          midSP.z;
          const auto rM =          midSP.r;
          const auto varianceZM =  midSP.varZ;
          const auto varianceRM =  midSP.varR;
          const auto cosPhiM =     xM / rM;
          const auto sinPhiM =     yM / rM;

          const auto deltaX = topSP.x - xM;
          const auto deltaY = topSP.y - yM;
          const auto deltaZ = topSP.z - zM;

          const auto x = deltaX * cosPhiM + deltaY * sinPhiM;
          const auto y = deltaY * cosPhiM - deltaX * sinPhiM;
          const auto iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
          const auto iDeltaR = cl::sycl::sqrt(iDeltaR2);
          const auto cot_theta = deltaZ * iDeltaR;

          deviceLinEqCircle L;
          L.cotTheta = cot_theta;
          L.zo = zM - rM * cot_theta;
          L.iDeltaR = iDeltaR;
          L.u = x * iDeltaR2;
          L.v = y * iDeltaR2;
          L.er = ((varianceZM + topSP.varZ) +
          (cot_theta * cot_theta) * (varianceRM + topSP.varR)) * iDeltaR2;

          linTopArray[idx] = L;
        });
      });

      //************************************************//
      // **** LINEAR EQUATION TRANSFORMATION - END **** //
      //************************************************//

      //************************************************//
      // *********** TRIPLET SEARCH - BEGIN *********** //
      //************************************************//

      const size_t maxMemoryAllocation = std::min(edgesComb, globalBufferSize / ((sizeof(TripletData)+sizeof(SeedData))*2));

      TripletData* curvImpactArray = cl::sycl::malloc_device<TripletData>(maxMemoryAllocation, (*q));
      SeedData * seedArray = cl::sycl::malloc_device<SeedData>(maxMemoryAllocation,*q);

      seeds.resize(M);
      const float MIN = -100000.f;
      size_t last_middle = 0;
      for(size_t first_middle = 0; first_middle < M; first_middle = last_middle){

        while(last_middle + 1 <= M &&
          (sumBotTopCombined[last_middle+1] - sumBotTopCombined[first_middle] < maxMemoryAllocation)){
            ++last_middle;
          }

        const uint32_t num_combinations = sumBotTopCombined[last_middle] - sumBotTopCombined[first_middle];
        if(num_combinations == 0) continue;

        const auto num_middle = last_middle - first_middle;

        uint32_t* offset = cl::sycl::malloc_device<uint32_t>(1,*q);
        std::atomic_uint* countTriplets = cl::sycl::malloc_device<std::atomic_uint>(1,*q); 
        std::atomic_uint* countSeeds = cl::sycl::malloc_device<std::atomic_uint>(1,*q); 
        q->memcpy(offset, &sumBotTopCombined[first_middle], sizeof(uint32_t));
        q->memset(countTriplets,0,sizeof(std::atomic_uint));
        q->memset(countSeeds,0,sizeof(std::atomic_uint));

        const auto sizeGlobalRange = (num_combinations/maxWorkItemPerGroup+1) *maxWorkItemPerGroup;
        cl::sycl::range<1> globalRange(sizeGlobalRange);
        cl::sycl::range<1> localRange(maxWorkItemPerGroup);        

        auto tripletEvent = q->submit([&](cl::sycl::handler &h) {
          h.depends_on({linB, linT});
          h.parallel_for<triplet_search_kernel>
            (cl::sycl::nd_range<1>{globalRange, localRange}, [=](cl::sycl::nd_item<1> item){
            const size_t idx = item.get_global_linear_id();
            if(idx < num_combinations){
            
              size_t L = first_middle;
              size_t R = last_middle;
              size_t mid = L;
              while(L < R - 1) {
                mid = (L + R) / 2;
                if(idx + offset[0] < sumCombArray[mid]) R = mid;
                else L = mid;
              }
              mid = L;

              TripletData T = {MIN, MIN};
              curvImpactArray[idx] = T;
              const auto numT = numTopArray[mid].fetch_add(0);

              const auto ib = sumBotArray[mid] + ((idx - sumCombArray[mid] + offset[0]) / numT);
              const auto it = sumTopArray[mid] + ((idx - sumCombArray[mid] + offset[0]) % numT);

              const auto bot = indBotArray[ib];
              const auto top = indTopArray[it];

              const auto linBotEq = linBotArray[ib];
              const auto linTopEq = linTopArray[it];
              const auto midSP =    midSPArray[mid];

              const auto Vb =          linBotEq.v;
              const auto Ub =          linBotEq.u;
              const auto Erb =         linBotEq.er;
              const auto cotThetab =   linBotEq.cotTheta;
              const auto iDeltaRb =    linBotEq.iDeltaR;

              const auto Vt =          linTopEq.v;
              const auto Ut =          linTopEq.u;
              const auto Ert =         linTopEq.er;
              const auto cotThetat =   linTopEq.cotTheta;
              const auto iDeltaRt =    linTopEq.iDeltaR;

              const auto rM =          midSP.r;
              const auto varianceRM =  midSP.varR;
              const auto varianceZM =  midSP.varZ;

              auto iSinTheta2 = (1. + cotThetab * cotThetab);
              auto scatteringInRegion2 = configData.maxScatteringAngle2 * iSinTheta2;
              scatteringInRegion2 *= configData.sigmaScattering * configData.sigmaScattering;
              auto error2 = Ert + Erb + 2 * (cotThetab * cotThetat *
                varianceRM + varianceZM) * iDeltaRb * iDeltaRt;
              auto deltaCotTheta = cotThetab - cotThetat;
              auto deltaCotTheta2 = deltaCotTheta * deltaCotTheta;

              deltaCotTheta = cl::sycl::abs(deltaCotTheta);
              auto error = cl::sycl::sqrt(error2);
              auto dCotThetaMinusError2 = deltaCotTheta2 + error2 - 2 * deltaCotTheta * error;
              auto dU = Ut - Ub;

              if((!(deltaCotTheta2 - error2 > 0) || !(dCotThetaMinusError2 > scatteringInRegion2))
                  && !(dU == 0.)) {
                auto A = (Vt - Vb) / dU;
                auto S2 = 1. + A * A;
                auto B = Vb - A * Ub;
                auto B2 = B * B;

                auto iHelixDiameter2 = B2 / S2;
                auto pT2scatter = 4 * iHelixDiameter2 * configData.pT2perRadius;
                auto p2scatter = pT2scatter * iSinTheta2;
                auto Im = cl::sycl::abs((A - B * rM) * rM);

                if(!(S2 < B2 * configData.minHelixDiameter2) && 
                    !((deltaCotTheta2 - error2 > 0) &&
                    (dCotThetaMinusError2 > p2scatter * configData.sigmaScattering * configData.sigmaScattering)) &&
                    !(Im > configData.impactMax)) {
                  auto c = countTriplets[0].fetch_add(1);
                  T.curvature = B / std::sqrt(S2);
                  T.impact = Im;
                  curvImpactArray[idx] = T;
                }
              }
            }
          });
        });
        tripletEvent.wait();

        uint sumTriplets;
        auto e0 = q->memcpy(&sumTriplets, countTriplets, sizeof(std::atomic_uint));
        e0.wait();

        if(sumTriplets == 0) continue;

        q->submit([&](cl::sycl::handler &h) {
          h.depends_on(tripletEvent);
          h.parallel_for<filter_2sp_fixed_kernel>
            (cl::sycl::nd_range<1>{globalRange, localRange}, [=](cl::sycl::nd_item<1> item){
            const size_t idx = item.get_global_linear_id();
            if(idx < num_combinations && curvImpactArray[idx].curvature != MIN){

              size_t L = first_middle;
              size_t R = last_middle;
              size_t mid = L;
              while(L < R - 1) {
                mid = (L + R) / 2;
                if(idx + offset[0] < sumCombArray[mid]) R = mid;
                else L = mid;
              }
              mid = L;

              const auto sumMidComb = sumCombArray[mid] - offset[0];
              const auto idxOffset = idx - sumMidComb;
              const auto numTopMid = numTopArray[mid].fetch_add(0);
              // const auto numTopMid = 20;

              const auto ib = sumBotArray[mid] + ((idxOffset) / numTopMid);
              const auto it = sumTopArray[mid] + ((idxOffset) % numTopMid);

              const auto bot = indBotArray[ib];
              const auto top = indTopArray[it];

              const auto current = curvImpactArray[idx];

              // by default compatSeedLimit is 2 -> 2 is hard coded
              // Variable length arrays are not supported in SYCL kernels.
              float compatibleSeedR[2];

              const auto invHelixDiameter = current.curvature;
              const auto lowerLimitCurv = invHelixDiameter - configData.deltaInvHelixDiameter;
              const auto upperLimitCurv = invHelixDiameter + configData.deltaInvHelixDiameter;
              const auto currentTop_r = topSPArray[top].r;
              auto weight = -(current.impact * configData.impactWeightFactor);

              size_t compatCounter = 0;

              const auto bottomOffset = ((idxOffset) / numTopMid) * numTopMid;

              for(uint32_t j = 0; j < numTopMid && compatCounter < configData.compatSeedLimit; ++j){
                uint32_t other_idx = sumMidComb + bottomOffset + j;
                float otherCurv = curvImpactArray[other_idx].curvature;
                if(otherCurv != MIN && other_idx != idx) {
                  uint32_t other_it = sumTopArray[mid] + j;
                  float otherTop_r = topSPArray[indTopArray[other_it]].r;
                  float deltaR = cl::sycl::abs(currentTop_r - otherTop_r);
                  if(deltaR >= configData.filterDeltaRMin &&
                    otherCurv >= lowerLimitCurv &&
                    otherCurv <= upperLimitCurv)
                    {
                      uint32_t c = 0;
                      for(; c < compatCounter &&
                          cl::sycl::abs(compatibleSeedR[c] - otherTop_r) >= configData.filterDeltaRMin; ++c){
                      }
                      if(c == compatCounter) {
                        compatibleSeedR[c] = otherTop_r;
                        ++compatCounter;
                      }
                    }
                }
              }

              weight += compatCounter * configData.compatSeedWeight;

              const auto bottomSP = botSPArray[bot];
              const auto middleSP = midSPArray[mid];
              const auto topSP = topSPArray[top];

              weight += deviceCuts.seedWeight(bottomSP, middleSP, topSP);

              if(deviceCuts.singleSeedCut(weight, bottomSP, middleSP, topSP)){
                const auto i = countSeeds[0].fetch_add(1);
                SeedData D;
                D.bottom = bot;
                D.top = top;
                D.middle = mid;
                D.weight = weight;
                seedArray[i] = D;
              }
            }
          });
        }).wait();

        uint sumSeeds;
        auto e1 = q->memcpy(&sumSeeds, countSeeds, sizeof(std::atomic_uint));
        e1.wait();
        
        std::vector<SeedData> hostSeedArray(sumSeeds);
        q->memcpy(&hostSeedArray[0], seedArray, sumSeeds* sizeof(SeedData));
        q->wait();

        for(uint32_t t = 0; t < sumSeeds; ++t) {
          auto m = hostSeedArray[t].middle;
          seeds[m].push_back(hostSeedArray[t]);
        }
      }

      //************************************************//
      // ************ TRIPLET SEARCH - END ************ //
      //************************************************//

      cl::sycl::free(tmpIndBotArray, *q);
      cl::sycl::free(tmpIndTopArray, *q);
      cl::sycl::free(botSPArray, *q);
      cl::sycl::free(midSPArray, *q);
      cl::sycl::free(topSPArray, *q);
      cl::sycl::free(numBotArray,*q);
      cl::sycl::free(numTopArray,*q);

      cl::sycl::free(linBotArray,*q);
      cl::sycl::free(linTopArray,*q);

      cl::sycl::free(indBotArray,*q);
      cl::sycl::free(indTopArray,*q);
      cl::sycl::free(indMidBotArray,*q);
      cl::sycl::free(indMidTopArray,*q);
      cl::sycl::free(sumBotArray,*q);
      cl::sycl::free(sumTopArray,*q);
      cl::sycl::free(sumCombArray, *q);
      
      cl::sycl::free(curvImpactArray,*q);
      cl::sycl::free(seedArray,*q);
    }
    catch (cl::sycl::exception const& e) {
      std::cout << "Caught synchronous SYCL exception:\n" << e.what() << std::endl;
      exit(0);
    }
  };
} // namespace Acts::Sycl
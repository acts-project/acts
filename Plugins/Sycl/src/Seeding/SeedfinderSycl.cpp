// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Sycl/Seeding/Seedfinder.hpp"
#include <CL/sycl.hpp>
#include <cstdint>
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

  using namespace Acts::Sycl::detail;

  void offloadComputations(cl::sycl::queue* q,
                          const deviceSeedfinderConfig& configData,
                          const DeviceExperimentCuts& deviceCuts,
                          const std::vector<deviceSpacePoint>& bottomSPs,
                          const std::vector<deviceSpacePoint>& middleSPs,
                          const std::vector<deviceSpacePoint>& topSPs,
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

    // Store the number of compatible bottom/top space points per middle space point.
    std::vector<uint32_t> numBotCompMid(M,0);
    std::vector<uint32_t> numTopCompMid(M,0);

    // Up to the Nth space point, the sum of compatible bottom/top space points.
    // We need these for indexing other vectors later in the algorithm.
    std::vector<uint32_t> sumBotCompUptoMid(M+1,0);
    std::vector<uint32_t> sumTopCompUptoMid(M+1,0);
    std::vector<uint64_t> sumBotTopCombined(M+1,0);

    // Similarly to indBotCompMid and indTopCompMid, we store the indices of the 
    // middle space points of the corresponding edges.
    std::vector<uint32_t> indMidBotComp;
    std::vector<uint32_t> indMidTopComp;

    // Number of edges for middle-bottom and middle-top duplet bipartite graphs.
    uint32_t edgesBottom = 0;
    uint32_t edgesTop = 0;
    uint64_t edgesComb = 0;

    // For the triplet search, we initialize some buffers outside the loop
    // for performance reasons. 
    uint32_t maxBotCompMid = 0;
    uint32_t maxTopCompMid = 0;

    try {
      using am = cl::sycl::access::mode;
      using at = cl::sycl::access::target;

      // Reserve buffers:
      //  - configBuf: device data of the SeedfinderConfig class instance m_config
      //  - botSPBuf, midSPBuf, topSPBuf: space point data
      //  - numBotCompBuf, numTopCompBuf: number of compatible bottom/top space points per middle sp
      cl::sycl::buffer<deviceSeedfinderConfig,1>     configBuf (&configData,1);
      cl::sycl::buffer<deviceSpacePoint,1>   botSPBuf  (bottomSPs.data(),  cl::sycl::range<1>(bottomSPs.size()));
      cl::sycl::buffer<deviceSpacePoint,1>   midSPBuf  (middleSPs.data(),  cl::sycl::range<1>(middleSPs.size()));
      cl::sycl::buffer<deviceSpacePoint,1>   topSPBuf  (topSPs.data(),     cl::sycl::range<1>(topSPs.size()));
      cl::sycl::buffer<uint32_t,1>     numBotCompBuf(numBotCompMid.data(), (cl::sycl::range<1>(M)));
      cl::sycl::buffer<uint32_t,1>     numTopCompBuf(numTopCompMid.data(), (cl::sycl::range<1>(M)));

      //*********************************************//
      // ********** DUPLET SEARCH - BEGIN ********** //
      //*********************************************//

      // The limit of compatible bottom [top] space points per middle space point is B [T].
      // Temporarily we reserve buffers of this size (M*B and M*T).
      // We store the indices of bottom [top] space points in bottomSPs [topSPs].
      // We move the indices to optimal size vectors for easier indexing.

      // After completing the duplet search, we'll have successfully contructed two
      // bipartite graphs for bottom-middle and top-middle space points.
      // We store the indices of the bottom/top space points of the edges of the graphs.
      // They index the bottomSPs and topSPs vectors.

      cl::sycl::buffer<uint32_t,1> tmpIndBotCompBuf((cl::sycl::range<1>(M*B)));
      cl::sycl::buffer<uint32_t,1> tmpIndTopCompBuf((cl::sycl::range<1>(M*T)));

      q->submit([&](cl::sycl::handler &h) {
        // Add accessors to buffers:
        auto configAcc =        configBuf.get_access<       am::read,           at::constant_buffer>(h);
        auto indBotCompatAcc =  tmpIndBotCompBuf.get_access<am::discard_write,  at::global_buffer>(h);
        auto numBotCompAcc =    numBotCompBuf.get_access<   am::atomic,         at::global_buffer>(h);
        auto botSPAcc =         botSPBuf.get_access<        am::read,           at::global_buffer>(h);
        auto midSPAcc =         midSPBuf.get_access<        am::read,           at::global_buffer>(h);

        h.parallel_for<duplet_search_bottom_kernel>
          (cl::sycl::range<2>{M,B}, [=](cl::sycl::id<2> idx) {
          const int mid = idx[0];
          const int bot = idx[1];

          deviceSpacePoint midSP = midSPAcc[mid];
          deviceSpacePoint botSP = botSPAcc[bot];
          deviceSeedfinderConfig config = configAcc[0];

          const float deltaR = midSP.r - botSP.r;
          const float cotTheta = (midSP.z - botSP.z) / deltaR;
          const float zOrigin = midSP.z - midSP.r * cotTheta;

          if( !(deltaR < config.deltaRMin) &&
              !(deltaR > config.deltaRMax) &&
              !(cl::sycl::abs(cotTheta) > config.cotThetaMax) &&
              !(zOrigin < config.collisionRegionMin) &&
              !(zOrigin > config.collisionRegionMax)
              && mid < M && bot < B) {
            // Besides checking the conditions that make a duplet compatible,
            // we also check whether this thread actually makes sense (within bounds).
            // We keep counting duplets with atomic access. 
            const int ind = numBotCompAcc[mid].fetch_add(1);
            indBotCompatAcc[mid * B + ind] = bot;
          }
        });
      });
  
      q->submit([&](cl::sycl::handler &h) {
        // Add accessors to buffers:
        auto configAcc =        configBuf.get_access<       am::read,           at::constant_buffer>(h);
        auto indTopCompatAcc =  tmpIndTopCompBuf.get_access<am::discard_write,  at::global_buffer>(h);
        auto numTopCompatAcc =  numTopCompBuf.get_access<   am::atomic,         at::global_buffer>(h);
        auto numBotCompAcc =    numBotCompBuf.get_access<   am::read,           at::global_buffer>(h);
        auto topSPAcc =         topSPBuf.get_access<        am::read,           at::global_buffer>(h);
        auto midSPAcc =         midSPBuf.get_access<        am::read,           at::global_buffer>(h);
        
        h.parallel_for<duplet_search_top_kernel>
          (cl::sycl::range<2>{M,B}, [=](cl::sycl::id<2> idx) {
          const int mid = idx[0];
          const int top = idx[1];

          deviceSpacePoint midSP = midSPAcc[mid];
          deviceSpacePoint topSP = topSPAcc[top];
          deviceSeedfinderConfig config = configAcc[0];

          if(numBotCompAcc[mid] != 0) {
            const float deltaR = topSP.r - midSP.r;
            const float cotTheta = (topSP.z - midSP.z) / deltaR;
            const float zOrigin = midSP.z - midSP.r * cotTheta;

            if( !(deltaR < config.deltaRMin) &&
                !(deltaR > config.deltaRMax) &&
                !(cl::sycl::abs(cotTheta) > config.cotThetaMax) &&
                !(zOrigin < config.collisionRegionMin) &&
                !(zOrigin > config.collisionRegionMax) &&
                mid < M && top < T) {
              // Besides checking the conditions that make a duplet compatible,
              // we also check whether this thread actually makes sense (within bounds).
              // We keep counting duplets with atomic access. 
              const int ind = numTopCompatAcc[mid].fetch_add(1);
              indTopCompatAcc[mid * T + ind] = top;
            }
          }
        });
      });

      //*********************************************//
      // *********** DUPLET SEARCH - END *********** //
      //*********************************************//

      // retrieve results from counting duplets
      {
        auto nB = numBotCompBuf.get_access<am::read>();
        auto nT = numTopCompBuf.get_access<am::read>();

        for(uint32_t i = 1; i < M + 1; ++i){
          sumBotCompUptoMid[i] += sumBotCompUptoMid[i-1] + nB[i-1];
          sumTopCompUptoMid[i] += sumTopCompUptoMid[i-1] + nT[i-1];
          sumBotTopCombined[i] += sumBotTopCombined[i-1] + nT[i-1]*nB[i-1];

          maxBotCompMid = std::max(maxBotCompMid, nB[i-1]);
          maxTopCompMid = std::max(maxTopCompMid, nT[i-1]);

          numBotCompMid[i-1] = nB[i-1];
          numTopCompMid[i-1] = nT[i-1]; 
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
      }

      // Copy indices from temporary matrices to final, optimal size vectors.
      // It will be simpler to index these buffers later on.

      cl::sycl::buffer<uint32_t,1> indTopCompBuf ((cl::sycl::range<1>(edgesTop)));
      cl::sycl::buffer<uint32_t,1> indBotCompBuf ((cl::sycl::range<1>(edgesBottom)));

      cl::sycl::buffer<uint32_t,1> indMidBotCompBuf (indMidBotComp.data(), cl::sycl::range<1>(edgesBottom));
      cl::sycl::buffer<uint32_t,1> indMidTopCompBuf (indMidTopComp.data(), cl::sycl::range<1>(edgesTop));

      cl::sycl::buffer<uint32_t,1> sumBotCompBuf (sumBotCompUptoMid.data(), cl::sycl::range<1>(M+1));
      cl::sycl::buffer<uint32_t,1> sumTopCompBuf (sumTopCompUptoMid.data(), cl::sycl::range<1>(M+1));      
      {
        q->submit([&](cl::sycl::handler &h){
          auto indBotAcc = indBotCompBuf.get_access<            am::write,  at::global_buffer>(h);
          auto indMidBotCompAcc = indMidBotCompBuf.get_access<  am::read,   at::global_buffer>(h);
          auto sumBotAcc = sumBotCompBuf.get_access<            am::read,   at::global_buffer>(h);
          auto tmpBotIndices = tmpIndBotCompBuf.get_access<     am::read,   at::global_buffer>(h);

          h.parallel_for<ind_copy_bottom_kernel>(edgesBottom, [=](cl::sycl::id<1> idx){
            uint32_t mid = indMidBotCompAcc[idx];
            uint32_t ind = tmpBotIndices[mid*B + idx - sumBotAcc[mid]];
            indBotAcc[idx] = ind;         
          });
        });

        q->submit([&](cl::sycl::handler &h){
          auto indTopAcc =          indTopCompBuf.get_access<     am::write,  at::global_buffer>(h);
          auto indMidTopCompAcc =   indMidTopCompBuf.get_access<  am::read,   at::global_buffer>(h);
          auto sumTopAcc =          sumTopCompBuf.get_access<     am::read,   at::global_buffer>(h);
          auto tmpTopIndices =      tmpIndTopCompBuf.get_access<  am::read,   at::global_buffer>(h);

          h.parallel_for<ind_copy_top_kernel>(edgesTop, [=](cl::sycl::id<1> idx){
            uint32_t mid = indMidTopCompAcc[idx];
            uint32_t ind = tmpTopIndices[mid*T + idx - sumTopAcc[mid]];
            
            indTopAcc[idx] = ind;     
          });
        });
      }    
    
      //************************************************//
      // *** LINEAR EQUATION TRANSFORMATION - BEGIN *** //
      //************************************************//

      // transformation of circle equation (x,y) into linear equation (u,v)
      // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
      // is transformed into
      // 1 - 2x_0*u - 2y_0*v = 0

      cl::sycl::buffer<deviceLinEqCircle,1> linBotBuf((cl::sycl::range<1>(edgesBottom)));
      cl::sycl::buffer<deviceLinEqCircle,1> linTopBuf((cl::sycl::range<1>(edgesTop)));

      // coordinate transformation middle-bottom pairs
      q->submit([&](cl::sycl::handler &h) {
        // add accessors to buffers
        auto indBotAcc =        indBotCompBuf.get_access<   am::read,         at::global_buffer>(h);
        auto indMidBotCompAcc = indMidBotCompBuf.get_access<am::read,         at::global_buffer>(h);
        auto sumBotAcc =        sumBotCompBuf.get_access<   am::read,         at::global_buffer>(h);
        auto botSPAcc =         botSPBuf.get_access<        am::read,         at::global_buffer>(h);
        auto midSPAcc =         midSPBuf.get_access<        am::read,         at::global_buffer>(h);
        auto linBotAcc =        linBotBuf.get_access<       am::discard_write,at::global_buffer>(h);

        h.parallel_for<transform_coord_bottom_kernel>(edgesBottom, [=](cl::sycl::id<1> idx) {
          deviceSpacePoint midSP = midSPAcc[indMidBotCompAcc[idx]];
          deviceSpacePoint botSP = botSPAcc[indBotAcc[idx]];

          float xM =          midSP.x;
          float yM =          midSP.y;
          float zM =          midSP.z;
          float rM =          midSP.r;
          float varianceZM =  midSP.varZ;
          float varianceRM =  midSP.varR;
          float cosPhiM =     xM / rM;
          float sinPhiM =     yM / rM;

          float deltaX = botSP.x - xM;
          float deltaY = botSP.y - yM;
          float deltaZ = botSP.z - zM;

          float x = deltaX * cosPhiM + deltaY * sinPhiM;
          float y = deltaY * cosPhiM - deltaX * sinPhiM;
          float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
          float iDeltaR = cl::sycl::sqrt(iDeltaR2);
          float cot_theta = -(deltaZ * iDeltaR);

          deviceLinEqCircle L;
          L.cotTheta = cot_theta;
          L.zo = zM - rM * cot_theta;
          L.iDeltaR = iDeltaR;
          L.u = x * iDeltaR2;
          L.v = y * iDeltaR2;
          L.er = ((varianceZM + botSP.varZ) +
          (cot_theta * cot_theta) * (varianceRM + botSP.varR)) * iDeltaR2;

          linBotAcc[idx] = L;
        });
      });
      
      // coordinate transformation middle-top pairs
      q->submit([&](cl::sycl::handler &h) {
        // add accessors to buffers
        auto indTopAcc =        indTopCompBuf.get_access<   am::read,         at::global_buffer>(h);
        auto indMidTopCompAcc = indMidTopCompBuf.get_access<am::read,         at::global_buffer>(h);
        auto sumTopAcc =        sumTopCompBuf.get_access<   am::read,         at::global_buffer>(h);
        auto topSPAcc =         topSPBuf.get_access<        am::read,         at::global_buffer>(h);
        auto midSPAcc =         midSPBuf.get_access<        am::read,         at::global_buffer>(h);
        auto linTopAcc =        linTopBuf.get_access<       am::discard_write,at::global_buffer>(h);

        h.parallel_for<transform_coord_top_kernel>(edgesTop, [=](cl::sycl::id<1> idx) {
          deviceSpacePoint midSP = midSPAcc[indMidTopCompAcc[idx]];
          deviceSpacePoint topSP = topSPAcc[indTopAcc[idx]];

          float xM =          midSP.x;
          float yM =          midSP.y;
          float zM =          midSP.z;
          float rM =          midSP.r;
          float varianceZM =  midSP.varZ;
          float varianceRM =  midSP.varR;
          float cosPhiM =     xM / rM;
          float sinPhiM =     yM / rM;

          float deltaX = topSP.x - xM;
          float deltaY = topSP.y - yM;
          float deltaZ = topSP.z - zM;

          float x = deltaX * cosPhiM + deltaY * sinPhiM;
          float y = deltaY * cosPhiM - deltaX * sinPhiM;
          float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
          float iDeltaR = cl::sycl::sqrt(iDeltaR2);
          float cot_theta = deltaZ * iDeltaR;

          deviceLinEqCircle L;
          L.cotTheta = cot_theta;
          L.zo = zM - rM * cot_theta;
          L.iDeltaR = iDeltaR;
          L.u = x * iDeltaR2;
          L.v = y * iDeltaR2;
          L.er = ((varianceZM + topSP.varZ) +
          (cot_theta * cot_theta) * (varianceRM + topSP.varR)) * iDeltaR2;

          linTopAcc[idx] = L;
        });
      });

      //************************************************//
      // **** LINEAR EQUATION TRANSFORMATION - END **** //
      //************************************************//

      //************************************************//
      // *********** TRIPLET SEARCH - BEGIN *********** //
      //************************************************//

      unsigned long globalBufferSize = q->get_device().get_info<cl::sycl::info::device::global_mem_size>();
      const uint32_t divide = 8;
      const size_t maxMemoryAllocation = std::min(edgesComb, globalBufferSize / (sizeof(TripletData)*divide));

      cl::sycl::buffer<uint64_t,1> sumCombinedBuf(sumBotTopCombined.data(), cl::sycl::range<1>(M+1));
      cl::sycl::buffer<TripletData,1> curvImpactBuf ((cl::sycl::range<1>(maxMemoryAllocation)));

      seeds.resize(M);
      const float MIN = -100000.f;
      uint32_t last_middle = 0;
      for(uint32_t first_middle = 0; first_middle < M; first_middle = last_middle){

        while(last_middle + 1 <= M &&
          (sumBotTopCombined[last_middle+1] - sumBotTopCombined[first_middle] < maxMemoryAllocation)){
            ++last_middle;
          }

        const uint64_t num_combinations = sumBotTopCombined[last_middle] - sumBotTopCombined[first_middle];

        if(num_combinations == 0) continue;

        const uint32_t num_bottoms = sumBotCompUptoMid[last_middle] - sumBotCompUptoMid[first_middle];
        const uint32_t num_tops = sumTopCompUptoMid[last_middle] - sumTopCompUptoMid[first_middle];
        const uint32_t num_middle = last_middle - first_middle;

        uint32_t countTriplets = 0;
        cl::sycl::buffer<uint32_t,1> countTripletsBuf(&countTriplets, 1);

        q->submit([&](cl::sycl::handler &h) {
          auto numTopAcc =      numTopCompBuf.get_access<   am::read,   at::global_buffer> (h);
          auto numBotAcc =      numBotCompBuf.get_access<   am::read,   at::global_buffer> (h);
          auto sumTopAcc =      sumTopCompBuf.get_access<   am::read,   at::global_buffer> (h);
          auto sumBotAcc =      sumBotCompBuf.get_access<   am::read,   at::global_buffer> (h);
          auto sumCombAcc =     sumCombinedBuf.get_access<  am::read,   at::global_buffer> (h);
          auto midSPAcc =       midSPBuf.get_access<        am::read,   at::global_buffer> (h);
          auto offsetComb =     sumCombinedBuf.get_access<  am::read,   at::constant_buffer> (h, 1, first_middle);
          auto offsetBot =      sumBotCompBuf.get_access<   am::read,   at::constant_buffer> (h, 1, first_middle);
          auto offsetTop =      sumTopCompBuf.get_access<   am::read,   at::constant_buffer> (h, 1, first_middle);
          auto configAcc =      configBuf.get_access<       am::read,   at::constant_buffer> (h);

          auto indBotAcc =      indBotCompBuf.get_access<   am::read,   at::global_buffer>
            (h, num_bottoms, sumBotCompUptoMid[first_middle]);
          auto indTopAcc =      indTopCompBuf.get_access<   am::read,   at::global_buffer>
            (h, num_tops, sumTopCompUptoMid[first_middle]);
          auto linBotAcc =      linBotBuf.get_access<       am::read,   at::global_buffer>
            (h, num_bottoms, sumBotCompUptoMid[first_middle]);
          auto linTopAcc =      linTopBuf.get_access<       am::read,   at::global_buffer>
            (h, num_tops, sumTopCompUptoMid[first_middle]);

          auto curvImpactAcc =  curvImpactBuf.get_access<   am::discard_write,at::global_buffer>
            (h, num_combinations, 0);
          auto countTripletsAcc=countTripletsBuf.get_access<am::atomic,       at::global_buffer>(h);

          h.parallel_for<triplet_search_kernel>
            (cl::sycl::range<1>{num_combinations}, [=](cl::sycl::id<1> idx){
            if(idx < num_combinations){
            
              uint32_t L = first_middle;
              uint32_t R = last_middle;
              uint32_t mid = L;
              while(L < R - 1) {
                mid = (L + R) / 2;
                if(idx + offsetComb[0] < sumCombAcc[mid]) R = mid;
                else L = mid;
              }
              mid = L;

              TripletData T = {MIN, MIN};
              curvImpactAcc[idx] = T;

              uint32_t ib = sumBotAcc[mid] - offsetBot[0] + ((idx - sumCombAcc[mid] + offsetComb[0]) / numTopAcc[mid]);
              uint32_t it = sumTopAcc[mid] - offsetTop[0] + ((idx - sumCombAcc[mid] + offsetComb[0]) % numTopAcc[mid]);

              uint32_t bot = indBotAcc[ib];
              uint32_t top = indTopAcc[it];
              deviceSeedfinderConfig config = configAcc[0];

              deviceLinEqCircle linBotEq = linBotAcc[ib];
              deviceLinEqCircle linTopEq = linTopAcc[it];
              deviceSpacePoint midSP = midSPAcc[mid];

              const float Vb =          linBotEq.v;
              const float Ub =          linBotEq.u;
              const float Erb =         linBotEq.er;
              const float cotThetab =   linBotEq.cotTheta;
              const float iDeltaRb =    linBotEq.iDeltaR;

              const float Vt =          linTopEq.v;
              const float Ut =          linTopEq.u;
              const float Ert =         linTopEq.er;
              const float cotThetat =   linTopEq.cotTheta;
              const float iDeltaRt =    linTopEq.iDeltaR;

              const float rM =          midSP.r;
              const float varianceRM =  midSP.varR;
              const float varianceZM =  midSP.varZ;

              float iSinTheta2 = (1. + cotThetab * cotThetab);
              float scatteringInRegion2 = config.maxScatteringAngle2 * iSinTheta2;
              scatteringInRegion2 *= config.sigmaScattering * config.sigmaScattering;
              float error2 = Ert + Erb + 2 * (cotThetab * cotThetat *
                varianceRM + varianceZM) * iDeltaRb * iDeltaRt;
              float deltaCotTheta = cotThetab - cotThetat;
              float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;

              deltaCotTheta = cl::sycl::abs(deltaCotTheta);
              float error = cl::sycl::sqrt(error2);
              float dCotThetaMinusError2 = deltaCotTheta2 + error2 - 2 * deltaCotTheta * error;

              float dU = Ut - Ub;

              if((!(deltaCotTheta2 - error2 > 0) || !(dCotThetaMinusError2 > scatteringInRegion2))
                  && !(dU == 0.)) {
                float A = (Vt - Vb) / dU;
                float S2 = 1. + A * A;
                float B = Vb - A * Ub;
                float B2 = B * B;

                float iHelixDiameter2 = B2 / S2;
                float pT2scatter = 4 * iHelixDiameter2 * config.pT2perRadius;
                float p2scatter = pT2scatter * iSinTheta2;
                float Im = cl::sycl::abs((A - B * rM) * rM);

                if(!(S2 < B2 * config.minHelixDiameter2) && 
                    !((deltaCotTheta2 - error2 > 0) &&
                    (dCotThetaMinusError2 > p2scatter * config.sigmaScattering * config.sigmaScattering)) &&
                    !(Im > config.impactMax)) {
                  uint32_t c = countTripletsAcc[0].fetch_add(1);
                  T.curvature = B / std::sqrt(S2);
                  T.impact = Im;
                  curvImpactAcc[idx] = T;
                }
              }
            }
          });
        });

        auto sumTriplets = (countTripletsBuf.get_access<am::read>())[0];
        if(sumTriplets == 0) continue;
        cl::sycl::buffer<SeedData,1> seedBuf((cl::sycl::range<1>(sumTriplets))); 

        uint32_t countSeeds = 0;
        cl::sycl::buffer<uint32_t,1>    countSeedsBuf(&countSeeds, 1);

        q->submit([&](cl::sycl::handler &h) {
          auto numTopAcc =      numTopCompBuf.get_access<   am::read,   at::global_buffer> (h);
          auto numBotAcc =      numBotCompBuf.get_access<   am::read,   at::global_buffer> (h);
          auto sumTopAcc =      sumTopCompBuf.get_access<   am::read,   at::global_buffer> (h);
          auto sumBotAcc =      sumBotCompBuf.get_access<   am::read,   at::global_buffer> (h);
          auto sumCombAcc =     sumCombinedBuf.get_access<  am::read,   at::global_buffer> (h);
          auto configAcc =      configBuf.get_access<       am::read,   at::constant_buffer> (h);
          auto offsetComb =     sumCombinedBuf.get_access<  am::read,   at::constant_buffer> (h, 1, first_middle);

          auto indBotAcc =      indBotCompBuf.get_access<   am::read,   at::global_buffer>(h);
          auto indTopAcc =      indTopCompBuf.get_access<   am::read,   at::global_buffer>(h);
          auto topSPAcc =       topSPBuf.get_access<        am::read,   at::global_buffer>(h);
          auto midSPAcc =       midSPBuf.get_access<        am::read,   at::global_buffer>(h);
          auto botSPAcc =       botSPBuf.get_access<        am::read,   at::global_buffer>(h);
          auto curvImpactAcc =  curvImpactBuf.get_access<   am::read,   at::global_buffer>(h, num_combinations, 0);

          auto seedAcc =        seedBuf.get_access<         am::write,  at::global_buffer>(h);
          auto countSeedsAcc=   countSeedsBuf.get_access<   am::atomic, at::global_buffer>(h);

          h.parallel_for<filter_2sp_fixed_kernel>
            (num_combinations, [=](cl::sycl::id<1> idx){
            if(idx < num_combinations && curvImpactAcc[idx].curvature != MIN){

              uint32_t L = first_middle;
              uint32_t R = last_middle;
              uint32_t mid = L;
              while(L < R - 1) {
                mid = (L + R) / 2;
                if(idx + offsetComb[0] < sumCombAcc[mid]) R = mid;
                else L = mid;
              }
              mid = L;

              uint32_t sumMidComb = sumCombAcc[mid] - offsetComb[0];
              uint32_t idxOffset = idx - sumMidComb;
              uint32_t numTopMid = numTopAcc[mid];

              uint32_t ib = sumBotAcc[mid] + ((idxOffset) / numTopMid);
              uint32_t it = sumTopAcc[mid] + ((idxOffset) % numTopMid);

              uint32_t bot = indBotAcc[ib];
              uint32_t top = indTopAcc[it];
              deviceSeedfinderConfig config = configAcc[0];

              TripletData current = curvImpactAcc[idx];

              // by default compatSeedLimit is 2 -> 2 is hard coded
              // Variable length arrays are not supported in SYCL kernels.
              float compatibleSeedR[2];

              float invHelixDiameter = current.curvature;
              float lowerLimitCurv = invHelixDiameter - config.deltaInvHelixDiameter;
              float upperLimitCurv = invHelixDiameter + config.deltaInvHelixDiameter;
              float currentTop_r = topSPAcc[top].r;
              float weight = -(current.impact * config.impactWeightFactor);

              uint32_t compatCounter = 0;

              uint32_t bottomOffset = ((idxOffset) / numTopMid) * numTopMid;

              for(uint32_t j = 0; j < numTopMid && compatCounter < config.compatSeedLimit; ++j){
                uint32_t other_idx = sumMidComb + bottomOffset + j;
                float otherCurv = curvImpactAcc[other_idx].curvature;
                if(otherCurv != MIN && other_idx != idx) {
                  uint32_t other_it = sumTopAcc[mid] + j;
                  float otherTop_r = topSPAcc[indTopAcc[other_it]].r;
                  float deltaR = cl::sycl::abs(currentTop_r - otherTop_r);
                  if(deltaR >= config.filterDeltaRMin &&
                    otherCurv >= lowerLimitCurv &&
                    otherCurv <= upperLimitCurv)
                    {
                      uint32_t c = 0;
                      for(; c < compatCounter &&
                          cl::sycl::abs(compatibleSeedR[c] - otherTop_r) >= config.filterDeltaRMin; ++c){
                      }
                      if(c == compatCounter) {
                        compatibleSeedR[c] = otherTop_r;
                        ++compatCounter;
                      }
                    }
                }
              }

              weight += compatCounter * config.compatSeedWeight;

              deviceSpacePoint bottomSP = botSPAcc[bot];
              deviceSpacePoint middleSP = midSPAcc[mid];
              deviceSpacePoint topSP = topSPAcc[top];

              weight += deviceCuts.seedWeight(bottomSP, middleSP, topSP);

              if(deviceCuts.singleSeedCut(weight, bottomSP, middleSP, topSP)){
                uint32_t i = countSeedsAcc[0].fetch_add(1);
                SeedData D;
                D.bottom = bot;
                D.top = top;
                D.middle = mid;
                D.weight = weight;
                seedAcc[i] = D;
              }
            }
          });
        });

        countSeeds = (countSeedsBuf.get_access<am::read>())[0];
        auto seedAcc = seedBuf.get_access<am::read>(countSeeds);

        for(uint32_t t = 0; t < countSeeds; ++t) {
          auto m = seedAcc[t].middle;
          seeds[m].push_back(seedAcc[t]);
        }
      }

      //************************************************//
      // ************ TRIPLET SEARCH - END ************ //
      //************************************************//
    }
    catch (cl::sycl::exception const& e) {
      std::cout << "Caught synchronous SYCL exception:\n" << e.what() << std::endl;
      exit(0);
    }
  };
} // namespace Acts::Sycl
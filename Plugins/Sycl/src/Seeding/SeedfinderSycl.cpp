// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Sycl/Seeding/Seedfinder.hpp"
#include <CL/sycl.hpp>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>
#include <array>
#include <exception>
#include <algorithm>


namespace Acts::Sycl {
  namespace sycl = cl::sycl;

  class triplet_search_kernel;
  class filter_2sp_fixed_kernel;

  void offloadComputations(cl::sycl::queue q,
                          const std::vector<float>& configData,
                          const std::vector<int>& limitData,
                          const std::vector<offloadSpacePoint>& bottomSPs,
                          const std::vector<offloadSpacePoint>& middleSPs,
                          const std::vector<offloadSpacePoint>& topSPs,
                          std::vector<std::vector<int>>& seedIndices,
                          std::vector<std::vector<float>>& seedWeight)
  {

    // Each vector stores data of space points in simplified 
    // structures of float variables
    // M: number of middle space points
    // B: number of bottom space points
    // T: number of top space points
    const size_t M = middleSPs.size();
    const size_t B = bottomSPs.size();
    const size_t T = topSPs.size();

    // Store the number of compatible bottom/top space points per middle space point.
    std::vector<int> numBotCompMid(M,0);
    std::vector<int> numTopCompMid(M,0);

    // Up to the Nth space point, the sum of compatible bottom/top space points.
    // We need these for indexing other vectors later in the algorithm.
    std::vector<int> sumBotCompUptoMid(M+1,0);
    std::vector<int> sumTopCompUptoMid(M+1,0);

    // After completing the duplet search, we'll have successfully contructed two
    // bipartite graphs for bottom-middle and top-middle space points.
    // We store the indices of the bottom/top space points of the edges of the graphs.
    // They index the bottomSPs and topSPs vectors.
    // std::vector<int> indBotCompMid;
    // std::vector<int> indTopCompMid;

    // Similarly to indBotCompMid and indTopCompMid, we store the indices of the 
    // middle space points of the corresponding edges.
    std::vector<int> indMidBotComp;
    std::vector<int> indMidTopComp;

    // Number of edges for middle-bottom and middle-top duplet bipartite graphs.
    int edgesBottom = 0;
    int edgesTop = 0;

    // For the triplet search, we initialize some buffers outside the loop
    // for performance reasons. 
    int maxBotCompMid = 0;
    int maxTopCompMid = 0;

    try {
      using am = sycl::access::mode;
      using at = sycl::access::target;

      // Reserve buffers:
      //  - configBuf: offloaded data of the SeedfinderConfig class instance m_config
      //  - limitBuf: required limit values, such as compatSeedLimit
      //  - botSPBuf, midSPBuf, topSPBuf: space point data
      //  - numBotCompBuf, numTopCompBuf: number of compatible bottom/top space points per middle sp
      sycl::buffer<float,1>   configBuf (configData.data(),         sycl::range<1>(configData.size()));
      sycl::buffer<int,1>     limitBuf  (limitData.data(),          sycl::range<1>(limitData.size()));
      sycl::buffer<offloadSpacePoint,1>
        botSPBuf  (bottomSPs.data(),          sycl::range<1>(bottomSPs.size()));
      sycl::buffer<offloadSpacePoint,1>
        midSPBuf  (middleSPs.data(),          sycl::range<1>(middleSPs.size()));
      sycl::buffer<offloadSpacePoint,1>
        topSPBuf  (topSPs.data(),             sycl::range<1>(topSPs.size()));
      sycl::buffer<int,1>     numBotCompBuf(numBotCompMid.data(), (sycl::range<1>(M)));
      sycl::buffer<int,1>     numTopCompBuf(numTopCompMid.data(), (sycl::range<1>(M)));

      //*********************************************//
      // ********** DUPLET SEARCH - BEGIN ********** //
      //*********************************************//

      // The limit of compatible bottom [top] space points per middle space point is B [T].
      // Temporarily we reserve buffers of this size (M*B and M*T).
      // We store the indices of bottom [top] space points in bottomSPs [topSPs].
      // We move the indices to optimal size vectors for algorithmic and performance reasons.

      sycl::buffer<int,1> tmpIndBotCompBuf((sycl::range<1>(M*B)));
      sycl::buffer<int,1> tmpIndTopCompBuf((sycl::range<1>(M*T)));

      // If it would be possible, we would only allocate memory for these buffers on the GPU side.
      // However, the current SYCL implementation reserves these on the host, and 

      auto bottom_duplet_search = q.submit([&](sycl::handler &h) {
        // Add accessors to buffers:
        auto configAcc =        configBuf.get_access<       am::read,           at::constant_buffer>(h);
        auto indBotCompatAcc =  tmpIndBotCompBuf.get_access<am::discard_write,  at::global_buffer>(h);
        auto numBotCompAcc =    numBotCompBuf.get_access<   am::atomic,         at::global_buffer>(h);
        auto botSPAcc =         botSPBuf.get_access<        am::read,           at::global_buffer>(h);
        auto midSPAcc =         midSPBuf.get_access<        am::read,           at::global_buffer>(h);

        h.parallel_for<class duplet_search_bottom>
          (sycl::range<2>{M,B}, [=](sycl::id<2> idx) {
          const int mid = idx[0];
          const int bot = idx[1];

          offloadSpacePoint midSP = midSPAcc[mid];
          offloadSpacePoint botSP = botSPAcc[bot];

          const float deltaR = midSP.r - botSP.r;
          const float cotTheta = (midSP.z - botSP.z) / deltaR;
          const float zOrigin = midSP.z - midSP.r * cotTheta;

          if( !(deltaR < configAcc[eDeltaRMin]) &&
              !(deltaR > configAcc[eDeltaRMax]) &&
              !(sycl::abs(cotTheta) > configAcc[eCotThetaMax]) &&
              !(zOrigin < configAcc[eCollisionRegionMin]) &&
              !(zOrigin > configAcc[eCollisionRegionMax])
              && mid < M && bot < B) {
            // Besides checking the conditions that make a duplet compatible,
            // we also check whether this thread actually makes sense (within bounds).
            // We keep counting duplets with atomic access. 
            const int ind = numBotCompAcc[mid].fetch_add(1);
            indBotCompatAcc[mid * B + ind] = bot;
          }
        });
      });
      bottom_duplet_search.wait();
  
      auto top_duplet_search = q.submit([&](sycl::handler &h) {
        // Add accessors to buffers:
        auto configAcc =        configBuf.get_access<       am::read,           at::constant_buffer>(h);
        auto indTopCompatAcc =  tmpIndTopCompBuf.get_access<am::discard_write,  at::global_buffer>(h);
        auto numTopCompatAcc =  numTopCompBuf.get_access<   am::atomic,         at::global_buffer>(h);
        auto numBotCompAcc =    numBotCompBuf.get_access<   am::read,           at::global_buffer>(h);
        auto topSPAcc =         topSPBuf.get_access<        am::read,           at::global_buffer>(h);
        auto midSPAcc =         midSPBuf.get_access<        am::read,           at::global_buffer>(h);
        
        h.parallel_for<class duplet_search_top>
          (sycl::range<2>{M,B}, [=](sycl::id<2> idx) {
          const int mid = idx[0];
          const int top = idx[1];

          offloadSpacePoint midSP = midSPAcc[mid];
          offloadSpacePoint topSP = topSPAcc[top];

          if(numBotCompAcc[mid] != 0) {
            const float deltaR = topSP.r - midSP.r;
            const float cotTheta = (topSP.z - midSP.z) / deltaR;
            const float zOrigin = midSP.z - midSP.r * cotTheta;

            if( !(deltaR < configAcc[eDeltaRMin]) &&
                !(deltaR > configAcc[eDeltaRMax]) &&
                !(sycl::abs(cotTheta) > configAcc[eCotThetaMax]) &&
                !(zOrigin < configAcc[eCollisionRegionMin]) &&
                !(zOrigin > configAcc[eCollisionRegionMax]) &&
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
      top_duplet_search.wait();

      //*********************************************//
      // *********** DUPLET SEARCH - END *********** //
      //*********************************************//

      // retrieve results from counting duplets
      {
        auto nB = numBotCompBuf.get_access<am::read>();
        auto nT = numTopCompBuf.get_access<am::read>();

        for(int i = 1; i < M + 1; ++i){
          sumBotCompUptoMid[i] += sumBotCompUptoMid[i-1] + nB[i-1];
          sumTopCompUptoMid[i] += sumTopCompUptoMid[i-1] + nT[i-1];

          maxBotCompMid = std::max(maxBotCompMid, nB[i-1]);
          maxTopCompMid = std::max(maxTopCompMid, nT[i-1]);

          numBotCompMid[i-1] = nB[i-1];
          numTopCompMid[i-1] = nT[i-1]; 
        }

        edgesBottom = sumBotCompUptoMid[M];
        edgesTop = sumTopCompUptoMid[M];

        if(edgesBottom == 0 || edgesTop == 0) return;

        indMidBotComp.reserve(edgesBottom);
        indMidBotComp.reserve(edgesTop);

        

        for(int mid = 0; mid < M; ++mid) {
          std::fill_n(std::back_inserter(indMidBotComp), nB[mid], mid);
          std::fill_n(std::back_inserter(indMidTopComp), nT[mid], mid);
        }
      }

      std::vector<int>    indTopCompMid(edgesTop);
      sycl::buffer<int,1> indBotCompBuf ((sycl::range<1>(edgesBottom)));

      sycl::buffer<int,1> indMidBotCompBuf (indMidBotComp.data(), sycl::range<1>(edgesBottom));
      sycl::buffer<int,1> indMidTopCompBuf (indMidTopComp.data(), sycl::range<1>(edgesTop));

      sycl::buffer<int,1> sumBotCompBuf (sumBotCompUptoMid.data(), sycl::range<1>(M+1));
      sycl::buffer<int,1> sumTopCompBuf (sumTopCompUptoMid.data(), sycl::range<1>(M+1));

      // Copy indices from temporary matrices to final, optimal size vectors.
      {
        sycl::buffer<int,1> indTopCompBuf (indTopCompMid.data(), sycl::range<1>(edgesTop));
        q.submit([&](sycl::handler &h){
          auto indBotAcc = indBotCompBuf.get_access<            am::write,  at::global_buffer>(h);
          auto indMidBotCompAcc = indMidBotCompBuf.get_access<  am::read,   at::global_buffer>(h);
          auto sumBotAcc = sumBotCompBuf.get_access<            am::read,   at::global_buffer>(h);
          auto tmpBotIndices = tmpIndBotCompBuf.get_access<     am::read,   at::global_buffer>(h);

          h.parallel_for<class ind_copy_bottom>(edgesBottom, [=](sycl::id<1> idx){
            int mid = indMidBotCompAcc[idx];
            int ind = tmpBotIndices[mid*B + idx - sumBotAcc[mid]];
            indBotAcc[idx] = ind;         
          });
        });

        q.submit([&](sycl::handler &h){
          auto indTopAcc =          indTopCompBuf.get_access<     am::write,  at::global_buffer>(h);
          auto indMidTopCompAcc =   indMidTopCompBuf.get_access<  am::read,   at::global_buffer>(h);
          auto sumTopAcc =          sumTopCompBuf.get_access<     am::read,   at::global_buffer>(h);
          auto tmpTopIndices =      tmpIndTopCompBuf.get_access<  am::read,   at::global_buffer>(h);

          h.parallel_for<class ind_copy_top>(edgesTop, [=](sycl::id<1> idx){
            int mid = indMidTopCompAcc[idx];
            int ind = tmpTopIndices[mid*T + idx - sumTopAcc[mid]];
            
            indTopAcc[idx] = ind;     
          });
        });
      } // indTopCompBuf buffer gets destroyed, data is copied back to indTopCompMid       

      // Sort by top indices for later filter algorithm. -> see filter_2sp_fixed_kernel
      // (ascending indices correspond to ascending radius because top space points are
      // already sorted by radius)
      for(int mid = 0; mid < M; ++mid){
        int sort_begin = sumTopCompUptoMid[mid];
        int sort_end = sumTopCompUptoMid[mid+1];
        if(sort_begin != sort_end) {
          std::sort(indTopCompMid.begin() + sort_begin, indTopCompMid.begin() + sort_end);
        }
      }
    
      //************************************************//
      // *** LINEAR EQUATION TRANSFORMATION - BEGIN *** //
      //************************************************//

      // transformation of circle equation (x,y) into linear equation (u,v)
      // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
      // is transformed into
      // 1 - 2x_0*u - 2y_0*v = 0

      sycl::buffer<int,1> indTopCompBuf (indTopCompMid.data(), sycl::range<1>(edgesTop));
      sycl::buffer<offloadLinEqCircle,1> linBotBuf((sycl::range<1>(edgesBottom)));
      sycl::buffer<offloadLinEqCircle,1> linTopBuf((sycl::range<1>(edgesTop)));

      // coordinate transformation middle-bottom pairs
      q.submit([&](sycl::handler &h) {
        // add accessors to buffers
        auto indBotAcc =        indBotCompBuf.get_access<   am::read,         at::global_buffer>(h);
        auto indMidBotCompAcc = indMidBotCompBuf.get_access<am::read,         at::global_buffer>(h);
        auto sumBotAcc =        sumBotCompBuf.get_access<   am::read,         at::global_buffer>(h);
        auto botSPAcc =         botSPBuf.get_access<        am::read,         at::global_buffer>(h);
        auto midSPAcc =         midSPBuf.get_access<        am::read,         at::global_buffer>(h);
        auto linBotAcc =        linBotBuf.get_access<       am::discard_write,at::global_buffer>(h);

        h.parallel_for<class transform_coord_bottom>(edgesBottom, [=](sycl::id<1> idx) {
          offloadSpacePoint midSP = midSPAcc[indMidBotCompAcc[idx]];
          offloadSpacePoint botSP = botSPAcc[indBotAcc[idx]];

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
          float iDeltaR = sycl::sqrt(iDeltaR2);
          float cot_theta = -(deltaZ * iDeltaR);

          offloadLinEqCircle L;
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
      q.submit([&](sycl::handler &h) {
        // add accessors to buffers
        auto indTopAcc =        indTopCompBuf.get_access<   am::read,         at::global_buffer>(h);
        auto indMidTopCompAcc = indMidTopCompBuf.get_access<am::read,         at::global_buffer>(h);
        auto sumTopAcc =        sumTopCompBuf.get_access<   am::read,         at::global_buffer>(h);
        auto topSPAcc =         topSPBuf.get_access<        am::read,         at::global_buffer>(h);
        auto midSPAcc =         midSPBuf.get_access<        am::read,         at::global_buffer>(h);
        auto linTopAcc =        linTopBuf.get_access<       am::discard_write,at::global_buffer>(h);

        h.parallel_for<class transform_coord_top>(edgesTop, [=](sycl::id<1> idx) {
          offloadSpacePoint midSP = midSPAcc[indMidTopCompAcc[idx]];
          offloadSpacePoint topSP = topSPAcc[indTopAcc[idx]];

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
          float iDeltaR = sycl::sqrt(iDeltaR2);
          float cot_theta = deltaZ * iDeltaR;

          offloadLinEqCircle L;
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

      seedIndices.resize(M);
      seedWeight.resize(M);
      const float MIN = -100000.f;
      sycl::buffer<float,3> curvImpactBuf ((sycl::range<3>(maxBotCompMid, maxTopCompMid, 2)));
      
      // Start kernels and load data to memory separately for each middle space point.
      // This way we don't run out of memory. (yay)
      for(int MID = 0; MID < M; ++MID) {
        if(numTopCompMid[MID] == 0 || numBotCompMid[MID] == 0) continue;

        // Count number of triplets for middle space point
        size_t maxTriplets = 0;
        sycl::buffer<size_t,1> maxTripletBuf(&maxTriplets, 1);

        q.submit([&](sycl::handler &h) {
          auto numTopAcc =      numTopCompBuf.get_access<   am::read,         at::global_buffer> (h, 1, MID);
          auto numBotAcc =      numBotCompBuf.get_access<   am::read,         at::global_buffer> (h, 1, MID);
          auto midSPAcc =       midSPBuf.get_access<        am::read,         at::global_buffer> (h, 1, MID);
          auto indBotAcc =      indBotCompBuf.get_access<   am::read,         at::global_buffer>
            (h, numBotCompMid[MID], sumBotCompUptoMid[MID]);
          auto indTopAcc =      indTopCompBuf.get_access<   am::read,         at::global_buffer>
            (h, numTopCompMid[MID], sumTopCompUptoMid[MID]);
          auto linBotAcc =      linBotBuf.get_access<       am::read,         at::global_buffer>
            (h, numBotCompMid[MID], sumBotCompUptoMid[MID]);
          auto linTopAcc =      linTopBuf.get_access<       am::read,         at::global_buffer>
            (h, numTopCompMid[MID], sumTopCompUptoMid[MID]);
          auto configAcc =      configBuf.get_access<       am::read,         at::constant_buffer>(h);

          auto curvImpactAcc =  curvImpactBuf.get_access<   am::discard_write,at::global_buffer>(h);
          auto maxTripletsAcc = maxTripletBuf.get_access<   am::atomic,       at::global_buffer>(h);

          h.parallel_for<triplet_search_kernel>
            (sycl::range<2>{size_t(numBotCompMid[MID]), size_t(numTopCompMid[MID])}, // number of threads
            [=](sycl::id<2> idx){
            // SYCL may start more threads than necessary (usually a power of 2)
            // which may cause us to find more triplets than what we should.
            // Check whether we are within bounds:
            // this costs us extra computing power, but gives better results.
            if(idx[0] < numBotAcc[0] && idx[1] < numTopAcc[0]) {
              curvImpactAcc[idx[0]][idx[1]][0] = MIN;
              curvImpactAcc[idx[0]][idx[1]][1] = MIN;

              int ib = idx[0];
              int it = idx[1];

              int bot = indBotAcc[ib];
              int top = indTopAcc[it];

              const float Zob =         linBotAcc[ib].zo;
              const float Vb =          linBotAcc[ib].v;
              const float Ub =          linBotAcc[ib].u;
              const float Erb =         linBotAcc[ib].er;
              const float cotThetab =   linBotAcc[ib].cotTheta;
              const float iDeltaRb =    linBotAcc[ib].iDeltaR;

              const float Zot =         linTopAcc[it].zo;
              const float Vt =          linTopAcc[it].v;
              const float Ut =          linTopAcc[it].u;
              const float Ert =         linTopAcc[it].er;
              const float cotThetat =   linTopAcc[it].cotTheta;
              const float iDeltaRt =    linTopAcc[it].iDeltaR;

              const float rM =          midSPAcc[0].r;
              const float varianceRM =  midSPAcc[0].varR;
              const float varianceZM =  midSPAcc[0].varZ;

              float iSinTheta2 = (1. + cotThetab * cotThetab);
              float scatteringInRegion2 = configAcc[eMaxScatteringAngle2] * iSinTheta2;
              scatteringInRegion2 *= configAcc[eSigmaScattering] * configAcc[eSigmaScattering];
              float error2 = Ert + Erb + 2 * (cotThetab * cotThetat *
                varianceRM + varianceZM) * iDeltaRb * iDeltaRt;
              float deltaCotTheta = cotThetab - cotThetat;
              float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;

              deltaCotTheta = sycl::abs(deltaCotTheta);
              float error = sycl::sqrt(error2);
              float dCotThetaMinusError2 = deltaCotTheta2 + error2 - 2 * deltaCotTheta * error;

              float dU = Ut - Ub;

              if((!(deltaCotTheta2 - error2 > 0) || !(dCotThetaMinusError2 > scatteringInRegion2))
                  && !(dU == 0.)) {
                float A = (Vt - Vb) / dU;
                float S2 = 1. + A * A;
                float B = Vb - A * Ub;
                float B2 = B * B;

                float iHelixDiameter2 = B2 / S2;
                float pT2scatter = 4 * iHelixDiameter2 * configAcc[ePT2perRadius];
                float p2scatter = pT2scatter * iSinTheta2;
                float Im = sycl::abs((A - B * rM) * rM);

                if(!(S2 < B2 * configAcc[eMinHelixDiameter2]) && 
                    !((deltaCotTheta2 - error2 > 0) &&
                    (dCotThetaMinusError2 > p2scatter * configAcc[eSigmaScattering] * configAcc[eSigmaScattering])) &&
                    !(Im > configAcc[eImpactMax])) {
                  maxTripletsAcc[0].fetch_add(1);
                  curvImpactAcc[idx[0]][idx[1]][0] = B / std::sqrt(S2);
                  curvImpactAcc[idx[0]][idx[1]][1] = Im;
                }
              }
            }
          });
        }).wait();
        maxTriplets = (maxTripletBuf.get_access<am::read>())[0];

        if(maxTriplets == 0) continue;

        // Reserve memory for triplet/seed indices and weights.
        // We could reserve these buffers outside the loop, and it would provide a performance
        // boost, if sub-buffers or sub-accessors had a correct implementation.
        // However, they don't.
        sycl::buffer<int,2>       seedIndBuf((sycl::range<2>{maxTriplets,2}));
        sycl::buffer<float,1>     seedWeightBuf((sycl::range<1>(maxTriplets)));        

        // Experiment specific cuts may reduce the number of triplets/seeds,
        // so we count them again. We only copy back to the host that many values.
        size_t countTriplets = 0;
        sycl::buffer<size_t,1>    countTripletsBuf(&countTriplets, 1);

        q.submit([&](sycl::handler &h) {
          // Since we only use part of the buffers, we can use sub-accessors.
          auto numTopAcc = numTopCompBuf.get_access<  am::read, at::global_buffer> (h, 1, MID);
          auto numBotAcc = numBotCompBuf.get_access<  am::read, at::global_buffer> (h, 1, MID);
          auto indBotAcc = indBotCompBuf.get_access<  am::read, at::global_buffer> (h, numBotCompMid[MID], sumBotCompUptoMid[MID]);
          auto indTopAcc = indTopCompBuf.get_access<  am::read, at::global_buffer> (h, numTopCompMid[MID], sumTopCompUptoMid[MID]);

          // Other buffers that we need to read data from
          auto curvImpactAcc =  curvImpactBuf.get_access<   am::read,         at::global_buffer>(h);
          auto topSPAcc =       topSPBuf.get_access<        am::read,         at::global_buffer>(h);
          auto botSPAcc =       botSPBuf.get_access<        am::read,         at::global_buffer>(h);
          auto configAcc =      configBuf.get_access<       am::read,         at::constant_buffer>(h);
          auto limitAcc =       limitBuf.get_access<        am::read,         at::constant_buffer>(h);

          // Use sub-accessors, so maybe with a future implementation this would
          // provide a benefit, currently it is no overhead
          auto seedIndAcc =  seedIndBuf.get_access<      am::write,        at::global_buffer>
            (h, sycl::range<2>{maxTriplets,2},sycl::range<2>{0,0});
          auto seedWeightAcc =  seedWeightBuf.get_access<   am::write,        at::global_buffer>
            (h,maxTriplets,0);
          auto countTripletsAcc=countTripletsBuf.get_access<am::atomic,       at::global_buffer>(h);

          h.parallel_for<filter_2sp_fixed_kernel>
            (sycl::range<2>{size_t(numBotCompMid[MID]), size_t(numTopCompMid[MID])}, // number of threads
            [=](sycl::id<2> idx){
            if(idx[0] < numBotAcc[0] && idx[1] < numTopAcc[0] 
                && curvImpactAcc[idx[0]][idx[1]][0] != MIN) {

              int bot = indBotAcc[idx[0]];
              int top = indTopAcc[idx[1]];

              float lowerLimitCurv = curvImpactAcc[idx[0]][idx[1]][0] - configAcc[eDeltaInvHelixDiameter];
              float upperLimitCurv = curvImpactAcc[idx[0]][idx[1]][0] + configAcc[eDeltaInvHelixDiameter];
              float currentTop_r = topSPAcc[top].r;
              float weight = -(curvImpactAcc[idx[0]][idx[1]][1] * configAcc[eImpactWeightFactor]);

              float lastCompatibleSeedR = currentTop_r;
              int compatCounter = 0;

              for(int j = 0; j < numTopAcc[0]; ++j){
                if(curvImpactAcc[idx[0]][j][0] != MIN && j != idx[1]) {
                  float otherTop_r = topSPAcc[indTopAcc[j]].r;
                  float deltaR = sycl::abs(currentTop_r - otherTop_r);
                  if(compatCounter < limitAcc[eCompatSeedLimit] &&
                    deltaR >= configAcc[eFilterDeltaRMin] &&
                    curvImpactAcc[idx[0]][j][0] >= lowerLimitCurv &&
                    curvImpactAcc[idx[0]][j][0] <= upperLimitCurv && 
                    sycl::abs(lastCompatibleSeedR - otherTop_r) >= configAcc[eFilterDeltaRMin]){
                      lastCompatibleSeedR = otherTop_r;
                      ++compatCounter;
                  }
                }
              }

              weight += compatCounter * configAcc[eCompatSeedWeight];

              // ATLAS experiment specific cuts
              float w = 0;
              if(botSPAcc[bot].r > 150){
                w = 400;
              }
              if(topSPAcc[top].r < 150){
                w = 200;
              }
              weight += w;

              if(!(botSPAcc[bot].r > 150. && weight < 380)) {
                int i = countTripletsAcc[0].fetch_add(1);
                seedIndAcc[i][0] = bot;
                seedIndAcc[i][1] = top;
                seedWeightAcc[i] = weight;
              }
            }
          });
        }).wait();
        countTriplets = (countTripletsBuf.get_access<am::read>())[0];

        // Sub-accessor to seed indices and weights.
        // Unfortunately this currently copies all data back
        // Hopefully, a future implementation would provide performance benefits
        auto seedInd = seedIndBuf.get_access<am::read>(sycl::range<2>{countTriplets,2},sycl::id<2>{0,0});
        auto seedWei = seedWeightBuf.get_access<am::read>(sycl::range<1>{countTriplets},sycl::id<1>{0});

        seedIndices[MID].reserve(countTriplets*2);
        seedWeight[MID].reserve(countTriplets);
        for(int j = 0; j < countTriplets; ++j) {
          seedIndices[MID].push_back(seedInd[j][0]);
          seedIndices[MID].push_back(seedInd[j][1]);
          seedWeight[MID].push_back(seedWei[j]);
        }
      }

      //************************************************//
      // ************ TRIPLET SEARCH - END ************ //
      //************************************************//
    }
    catch (sycl::exception const& e) {
      std::cout << "Caught synchronous SYCL exception:\n" << e.what() << std::endl;
    }
  };
} // namespace Acts::Sycl
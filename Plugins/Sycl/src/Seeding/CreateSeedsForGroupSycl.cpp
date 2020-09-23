// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// System include(s)
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <exception>
#include <functional>
#include <vector>

// Acts include(s)
#include "Acts/Utilities/Logger.hpp"

// SYCL plugin include(s)
#include "Acts/Plugins/Sycl/Seeding/CreateSeedsForGroupSycl.hpp"
#include "Acts/Plugins/Sycl/Seeding/detail/Types.hpp"
#include "Acts/Plugins/Sycl/Utilities/CalculateNdRange.hpp"

// SYCL include
#include <CL/sycl.hpp>

namespace Acts::Sycl {
// Kernel classes in order of execution.
class duplet_search_bottom_kernel;
class duplet_search_top_kernel;
class ind_copy_bottom_kernel;
class ind_copy_top_kernel;
class transform_coord_bottom_kernel;
class transform_coord_top_kernel;
class triplet_search_kernel;
class filter_2sp_fixed_kernel;

void createSeedsForGroupSycl(
    const QueueWrapper& wrappedQueue,
    const detail::DeviceSeedfinderConfig& seedfinderConfig,
    const DeviceExperimentCuts& deviceCuts,
    const std::vector<detail::DeviceSpacePoint>& bottomSPs,
    const std::vector<detail::DeviceSpacePoint>& middleSPs,
    const std::vector<detail::DeviceSpacePoint>& topSPs,
    std::vector<std::vector<detail::SeedData>>& seeds) {
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
  // These are prefix sum arrays, with a leading zero.
  std::vector<uint32_t> sumBotCompUptoMid(M + 1, 0);
  std::vector<uint32_t> sumTopCompUptoMid(M + 1, 0);
  std::vector<uint32_t> sumBotTopCombined(M + 1, 0);

  // After completing the duplet search, we'll have successfully contructed
  // two bipartite graphs for bottom-middle and top-middle space points.
  // We store the indices of the middle space points of the corresponding
  // edges.
  std::vector<uint32_t> indMidBotComp;
  std::vector<uint32_t> indMidTopComp;

  // Number of edges for middle-bottom and middle-top duplet bipartite graphs.
  uint64_t edgesBottom = 0;
  uint64_t edgesTop = 0;
  // Number of possible compatible triplets. This is the sum of the
  // combination of the number of compatible bottom and compatible top duplets
  // per middle space point. (nb0*nt0 + nb1*nt1 + ... where nbk is the number
  // of comp. bot. SPs for the kth middle SP)
  uint64_t edgesComb = 0;

  try {
    auto* q = wrappedQueue.getQueue();
    uint64_t globalBufferSize =
        q->get_device().get_info<cl::sycl::info::device::global_mem_size>();
    uint64_t maxWorkGroupSize =
        q->get_device().get_info<cl::sycl::info::device::max_work_group_size>();

    // Device allocations
    detail::DeviceSpacePoint* deviceBottomSPs =
        cl::sycl::malloc_device<detail::DeviceSpacePoint>(B, *q);
    detail::DeviceSpacePoint* deviceMiddleSPs =
        cl::sycl::malloc_device<detail::DeviceSpacePoint>(M, *q);
    detail::DeviceSpacePoint* deviceTopSPs =
        cl::sycl::malloc_device<detail::DeviceSpacePoint>(T, *q);

    // The limit of compatible bottom [top] space points per middle space
    // point is B [T]. Temporarily we reserve buffers of this size (M*B and
    // M*T). We store the indices of bottom [top] space points in bottomSPs
    // [topSPs]. We move the indices to optimal size vectors for easier
    // indexing.
    uint32_t* deviceTmpIndBot = cl::sycl::malloc_device<uint32_t>(M * B, *q);
    uint32_t* deviceTmpIndTop = cl::sycl::malloc_device<uint32_t>(M * T, *q);

    q->memcpy(deviceBottomSPs, bottomSPs.data(),
              sizeof(detail::DeviceSpacePoint) * (B));
    q->memcpy(deviceMiddleSPs, middleSPs.data(),
              sizeof(detail::DeviceSpacePoint) * (M));
    q->memcpy(deviceTopSPs, topSPs.data(),
              sizeof(detail::DeviceSpacePoint) * (T));

    // Calculate 2 dimensional range of bottom-middle duplet search kernel
    // We'll have a total of M*B threads globally, but we need to give the
    // nd_range the global dimensions so that they are an exact multiple of
    // the local dimensions. That's why we need this calculation.

    cl::sycl::nd_range<2> bottomDupletNDRange =
        calculate2DimNDRange(M, B, maxWorkGroupSize);
    cl::sycl::nd_range<2> topDupletNDRange =
        calculate2DimNDRange(M, T, maxWorkGroupSize);

    // Count compatible bottom/top space points per middle space point.
    std::vector<uint32_t> countBotDuplets(M, 0);
    std::vector<uint32_t> countTopDuplets(M, 0);

    //*********************************************//
    // ********** DUPLET SEARCH - BEGIN ********** //
    //*********************************************//

    {
      sycl::buffer<uint32_t> countBotBuf(countBotDuplets.data(), M);
      sycl::buffer<uint32_t> countTopBuf(countTopDuplets.data(), M);
      q->submit([&](cl::sycl::handler& h) {
        auto countBotDupletsAcc =
            sycl::ONEAPI::atomic_accessor<uint32_t, 1,
                                          sycl::ONEAPI::memory_order::relaxed,
                                          sycl::ONEAPI::memory_scope::device>(
                countBotBuf, h);
        h.parallel_for<duplet_search_bottom_kernel>(
            bottomDupletNDRange, [=](cl::sycl::nd_item<2> item) {
              const auto mid = item.get_global_id(0);
              const auto bot = item.get_global_id(1);
              // We check whether this thread actually makes sense (within
              // bounds).
              // The number of threads is usually a factor of 2, or 3*2^k (k
              // \in N), etc. Without this check we may index out of arrays.
              if (mid < M && bot < B) {
                const auto midSP = deviceMiddleSPs[mid];
                const auto botSP = deviceBottomSPs[bot];

                const auto deltaR = midSP.r - botSP.r;
                const auto cotTheta = (midSP.z - botSP.z) / deltaR;
                const auto zOrigin = midSP.z - midSP.r * cotTheta;

                if (!(deltaR < seedfinderConfig.deltaRMin) &&
                    !(deltaR > seedfinderConfig.deltaRMax) &&
                    !(cl::sycl::abs(cotTheta) > seedfinderConfig.cotThetaMax) &&
                    !(zOrigin < seedfinderConfig.collisionRegionMin) &&
                    !(zOrigin > seedfinderConfig.collisionRegionMax)) {
                  // We keep counting duplets with atomic access.
                  const auto ind = countBotDupletsAcc[mid].fetch_add(1);
                  deviceTmpIndBot[mid * B + ind] = bot;
                }
              }
            });
      });

      q->submit([&](cl::sycl::handler& h) {
        auto countTopDupletsAcc =
            sycl::ONEAPI::atomic_accessor<uint32_t, 1,
                                          sycl::ONEAPI::memory_order::relaxed,
                                          sycl::ONEAPI::memory_scope::device>(
                countTopBuf, h);
        h.parallel_for<duplet_search_top_kernel>(
            topDupletNDRange, [=](cl::sycl::nd_item<2> item) {
              const auto mid = item.get_global_id(0);
              const auto top = item.get_global_id(1);
              // We check whether this thread actually makes sense (within
              // bounds).
              if (mid < M && top < T) {
                const auto midSP = deviceMiddleSPs[mid];
                const auto topSP = deviceTopSPs[top];

                const auto deltaR = topSP.r - midSP.r;
                const auto cotTheta = (topSP.z - midSP.z) / deltaR;
                const auto zOrigin = midSP.z - midSP.r * cotTheta;

                if (!(deltaR < seedfinderConfig.deltaRMin) &&
                    !(deltaR > seedfinderConfig.deltaRMax) &&
                    !(cl::sycl::abs(cotTheta) > seedfinderConfig.cotThetaMax) &&
                    !(zOrigin < seedfinderConfig.collisionRegionMin) &&
                    !(zOrigin > seedfinderConfig.collisionRegionMax)) {
                  // We keep counting duplets with atomic access.
                  const auto ind = countTopDupletsAcc[mid].fetch_add(1);
                  deviceTmpIndTop[mid * T + ind] = top;
                }
              }
            });
      });
    }  // sync (buffers get destroyed and duplet counts are copied back to
       // countBotDuplets and countTopDuplets)

    //*********************************************//
    // *********** DUPLET SEARCH - END *********** //
    //*********************************************//

    // retrieve results from counting duplets
    {
      // Construct prefix sum arrays of duplet counts.
      // These will later be used to index other arrays based on middle SP
      // indices.
      for (uint32_t i = 1; i < M + 1; ++i) {
        sumBotCompUptoMid[i] +=
            sumBotCompUptoMid[i - 1] + countBotDuplets[i - 1];
        sumTopCompUptoMid[i] +=
            sumTopCompUptoMid[i - 1] + countTopDuplets[i - 1];
        sumBotTopCombined[i] += sumBotTopCombined[i - 1] +
                                countTopDuplets[i - 1] * countBotDuplets[i - 1];
      }

      edgesBottom = sumBotCompUptoMid[M];
      edgesTop = sumTopCompUptoMid[M];
      edgesComb = sumBotTopCombined[M];

      indMidBotComp.reserve(edgesBottom);
      indMidTopComp.reserve(edgesTop);

      // Fill arrays of middle SP indices of found duplets (bottom and top).
      for (uint32_t mid = 0; mid < M; ++mid) {
        std::fill_n(std::back_inserter(indMidBotComp), countBotDuplets[mid],
                    mid);
        std::fill_n(std::back_inserter(indMidTopComp), countTopDuplets[mid],
                    mid);
      }
    }

    if (edgesBottom > 0 && edgesTop > 0) {
      // Calcualte global and local range of execution for edgesBottom number of
      // threads. Local range is the same as block size in CUDA.
      cl::sycl::nd_range<1> edgesBotNdRange =
          calculate1DimNDRange(edgesBottom, maxWorkGroupSize);

      // Global and local range of execution for edgesTop number of threads.
      cl::sycl::nd_range<1> edgesTopNdRange =
          calculate1DimNDRange(edgesTop, maxWorkGroupSize);

      // EXPLANATION OF INDEXING (fisrt part)
      /*
        (for bottom-middle duplets, but it is the same for middle-tops)

        In case we have 4 middle SP and 5 bottom SP, our temporary array of
        the compatible bottom duplet indices would look like this:
             ---------------------
        mid0 | 0 | 3 | 4 | 1 | - |    Indices in the columns correspond to
        mid1 | 3 | 2 | - | - | - |    bottom SP indices in the bottomSPs
        mid2 | - | - | - | - | - |    array. Threads are executed concurrently,
        mid3 | 4 | 2 | 1 | - | - |    so the order of indices is random.
             ---------------------
        We will refer to this structure as a bipartite graph, as it can be
        described by a graph of nodes for middle and bottom SPs, and edges
        between one middle and one bottom SP, but never two middle or two
        bottom SPs.

        We will flatten this matrix out, and store the indices the
        following way (this is deviceIndBot):
        -------------------------------------
        | 0 | 3 | 4 | 1 | 3 | 2 | 4 | 2 | 1 |
        -------------------------------------

        Also the length of this array is equal to edgesBottom, which is 9 in
        this example. It is the number of the edges of the bottom-middle
        bipartite graph.

        To find out where the indices of bottom SPs start for a particular
        middle SP, we use prefix sum arrays.
        We know how many duplets were found for each middle SP (this is
        deviceCountBotDuplets).
        -----------------
        | 4 | 2 | 0 | 3 |
        -----------------

        We will make a prefix sum array of these counts, with a leading zero:
        (this is deviceSumBot)
        ---------------------
        | 0 | 4 | 6 | 6 | 9 |
        ---------------------

        If we have the middle SP with index 1, then we know that the indices
        of the compatible bottom SPs are in the range (left closed, right
        open) [deviceSumBot[1] , deviceSumBot[2] ) of deviceIndBot. In this
        case, these indices are 3 and 2, so we'd use these to index
        deviceBottomSPs to gather data about the bottom SP.

        To be able to get the indices of middle SPs in constant time inside
        kernels, we will also prepare arrays that store the indices of the
        middleSPs of the edges (deviceMidIndPerBot).
        -------------------------------------
        | 0 | 0 | 0 | 0 | 1 | 1 | 3 | 3 | 3 |
        -------------------------------------

        (For the same purpose, we could also do a binary search on the
        deviceSumBot array, and we will do exactly that later, in the triplet
        search kernel.)

        We will execute the coordinate transformation on edgesBottom threads,
        or 9 in our example.

        The size of the array storing our transformed coordinates
        (deviceLinBot) is also edgesBottom, the sum of bottom duplets we
        found so far.
      */

      sycl::buffer<uint32_t> numTopDupletsBuf(countTopDuplets.data(), M);

      // We store the indices of the BOTTOM/TOP space points of the edges of
      // the bottom-middle and top-middle bipartite duplet graphs. They index
      // the bottomSPs and topSPs vectors.
      uint32_t* deviceIndBot =
          cl::sycl::malloc_device<uint32_t>(edgesBottom, *q);
      uint32_t* deviceIndTop = cl::sycl::malloc_device<uint32_t>(edgesTop, *q);

      // We store the indices of the MIDDLE space points of the edges of the
      // bottom-middle and top-middle bipartite duplet graphs.
      // They index the middleSP vector.
      uint32_t* deviceMidIndPerBot =
          cl::sycl::malloc_device<uint32_t>(edgesBottom, *q);
      uint32_t* deviceMidIndPerTop =
          cl::sycl::malloc_device<uint32_t>(edgesTop, *q);

      // Partial sum arrays of deviceNumBot and deviceNum
      uint32_t* deviceSumBot = cl::sycl::malloc_device<uint32_t>(M + 1, *q);
      uint32_t* deviceSumTop = cl::sycl::malloc_device<uint32_t>(M + 1, *q);

      // Partial sum array of the combinations of compatible bottom and top
      // space points per middle space point.
      uint32_t* deviceSumComb = cl::sycl::malloc_device<uint32_t>(M + 1, *q);

      // Allocations for coordinate transformation.
      detail::DeviceLinEqCircle* deviceLinBot =
          cl::sycl::malloc_device<detail::DeviceLinEqCircle>(edgesBottom, *q);
      detail::DeviceLinEqCircle* deviceLinTop =
          cl::sycl::malloc_device<detail::DeviceLinEqCircle>(edgesTop, *q);

      q->memcpy(deviceMidIndPerBot, indMidBotComp.data(),
                sizeof(uint32_t) * edgesBottom);
      q->memcpy(deviceMidIndPerTop, indMidTopComp.data(),
                sizeof(uint32_t) * edgesTop);
      q->memcpy(deviceSumBot, sumBotCompUptoMid.data(),
                sizeof(uint32_t) * (M + 1));
      q->memcpy(deviceSumTop, sumTopCompUptoMid.data(),
                sizeof(uint32_t) * (M + 1));
      q->memcpy(deviceSumComb, sumBotTopCombined.data(),
                sizeof(uint32_t) * (M + 1));
      q->wait();

      // Copy indices from temporary matrices to final, optimal size vectors.
      // We will use these for easier indexing.
      {
        q->submit([&](cl::sycl::handler& h) {
          h.parallel_for<ind_copy_bottom_kernel>(
              edgesBotNdRange, [=](cl::sycl::nd_item<1> item) {
                auto idx = item.get_global_linear_id();
                if (idx < edgesBottom) {
                  auto mid = deviceMidIndPerBot[idx];
                  auto ind = deviceTmpIndBot[mid * B + idx - deviceSumBot[mid]];
                  deviceIndBot[idx] = ind;
                }
              });
        });

        q->submit([&](cl::sycl::handler& h) {
          h.parallel_for<ind_copy_top_kernel>(
              edgesTopNdRange, [=](cl::sycl::nd_item<1> item) {
                auto idx = item.get_global_linear_id();
                if (idx < edgesTop) {
                  auto mid = deviceMidIndPerTop[idx];
                  auto ind = deviceTmpIndTop[mid * T + idx - deviceSumTop[mid]];
                  deviceIndTop[idx] = ind;
                }
              });
        });
      }  // sync

      //************************************************//
      // *** LINEAR EQUATION TRANSFORMATION - BEGIN *** //
      //************************************************//

      // transformation of circle equation (x,y) into linear equation (u,v)
      // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
      // is transformed into
      // 1 - 2x_0*u - 2y_0*v = 0

      // coordinate transformation middle-bottom pairs
      auto linB = q->submit([&](cl::sycl::handler& h) {
        h.parallel_for<transform_coord_bottom_kernel>(
            edgesBotNdRange, [=](cl::sycl::nd_item<1> item) {
              auto idx = item.get_global_linear_id();
              if (idx < edgesBottom) {
                const auto midSP = deviceMiddleSPs[deviceMidIndPerBot[idx]];
                const auto botSP = deviceBottomSPs[deviceIndBot[idx]];

                const auto xM = midSP.x;
                const auto yM = midSP.y;
                const auto zM = midSP.z;
                const auto rM = midSP.r;
                const auto varianceZM = midSP.varZ;
                const auto varianceRM = midSP.varR;
                const auto cosPhiM = xM / rM;
                const auto sinPhiM = yM / rM;

                const auto deltaX = botSP.x - xM;
                const auto deltaY = botSP.y - yM;
                const auto deltaZ = botSP.z - zM;

                const auto x = deltaX * cosPhiM + deltaY * sinPhiM;
                const auto y = deltaY * cosPhiM - deltaX * sinPhiM;
                const auto iDeltaR2 = 1.f / (deltaX * deltaX + deltaY * deltaY);
                const auto iDeltaR = cl::sycl::sqrt(iDeltaR2);
                const auto cot_theta = -(deltaZ * iDeltaR);

                detail::DeviceLinEqCircle L;
                L.cotTheta = cot_theta;
                L.zo = zM - rM * cot_theta;
                L.iDeltaR = iDeltaR;
                L.u = x * iDeltaR2;
                L.v = y * iDeltaR2;
                L.er = ((varianceZM + botSP.varZ) +
                        (cot_theta * cot_theta) * (varianceRM + botSP.varR)) *
                       iDeltaR2;

                deviceLinBot[idx] = L;
              }
            });
      });

      // coordinate transformation middle-top pairs
      auto linT = q->submit([&](cl::sycl::handler& h) {
        h.parallel_for<transform_coord_top_kernel>(
            edgesTopNdRange, [=](cl::sycl::nd_item<1> item) {
              auto idx = item.get_global_linear_id();
              if (idx < edgesTop) {
                const auto midSP = deviceMiddleSPs[deviceMidIndPerTop[idx]];
                const auto topSP = deviceTopSPs[deviceIndTop[idx]];

                const auto xM = midSP.x;
                const auto yM = midSP.y;
                const auto zM = midSP.z;
                const auto rM = midSP.r;
                const auto varianceZM = midSP.varZ;
                const auto varianceRM = midSP.varR;
                const auto cosPhiM = xM / rM;
                const auto sinPhiM = yM / rM;

                const auto deltaX = topSP.x - xM;
                const auto deltaY = topSP.y - yM;
                const auto deltaZ = topSP.z - zM;

                const auto x = deltaX * cosPhiM + deltaY * sinPhiM;
                const auto y = deltaY * cosPhiM - deltaX * sinPhiM;
                const auto iDeltaR2 = 1.f / (deltaX * deltaX + deltaY * deltaY);
                const auto iDeltaR = cl::sycl::sqrt(iDeltaR2);
                const auto cot_theta = deltaZ * iDeltaR;

                detail::DeviceLinEqCircle L;
                L.cotTheta = cot_theta;
                L.zo = zM - rM * cot_theta;
                L.iDeltaR = iDeltaR;
                L.u = x * iDeltaR2;
                L.v = y * iDeltaR2;
                L.er = ((varianceZM + topSP.varZ) +
                        (cot_theta * cot_theta) * (varianceRM + topSP.varR)) *
                       iDeltaR2;

                deviceLinTop[idx] = L;
              }
            });
      });

      //************************************************//
      // **** LINEAR EQUATION TRANSFORMATION - END **** //
      //************************************************//

      //************************************************//
      // *********** TRIPLET SEARCH - BEGIN *********** //
      //************************************************//

      // EXPLANATION OF INDEXING (second part)
      /*
        For the triplet search, we calculate the upper limit of constructible
        triplets.

        For this, we multiply the number of compatible bottom and compatible
        top SPs for each middle SP, and add these together.
        (nb0*nt0 + nb1*nt1 + ... where nbk is the number of compatible bottom
        SPs for the kth middle SP, similarly ntb is for tops)

        sumBotTopCombined is a prefix sum array (of length M+1) of the
        calculated combinations.

        sumBotTopCombined:
        ________________________________________________________
        |     |         |                   |     |  M         | M = number
        |  0  | nb0*nt0 | nb0*nt0 + nb1*nt1 | ... |  ∑ nbi+nti | of middle
        |_____|_________|___________________|_____|_i=0________| space points

        We will start kernels and reserve memory for these combinations but
        only so much we can fit into memory at once.

        We limit our memory usage to globalBufferSize/2, this is currently
        hard-coded, but it could be configured. Actually, it would be better
        to use a separate object that manages memory allocations and
        deallocations and we could ask it to lend us as much memory as it is
        happy to give.

        For later, let maxMemoryAllocation be maximum allocatable memory for
        triplet search.

        We start by adding up summing the combinations, until we arrive at a
        k which for:

        k+1
         ∑ nbi+nti > maxMemoryAllocation
        i=0
        (or k == M).

        So we know, that we need to start our first kernel for the first k
        middle SPs.

        Inside the triplet search kernel we start with a binary search, to
        find out which middle SP the thread corresponds to. Note, that
        sumBotTopCombined is a monotone increasing series of values which
        allows us to do a binary search on it.

        Inside the triplet search kernel we count the triplets for fixed
        bottom and middle SP. This is deviceCountTriplets.

        The triplet filter kernel is calculated on threads equal to all possible
        bottom-middle combinations for the first k middle SPs, which are
        the sum of bottom-middle duplets. (For the next kernel it would be the
        bottom-middle combinations from the (k+1)th middle SP to another jth
        middle SP j<=M.)

        This will be numTripletFilterThreads =
            sumBotCompUptoMid[lastMiddle] - sumBotCompUptoMid[firstMiddle]

        If the triplet search and triplet filter kernel finished, we continue
        summing up possible triplet combinations from the (k+1)th middle SP.

        Inside the kernels we need to use offset because of this, to be able to
        map threads to space point indices.

        This offset is sumCombUptoFirstMiddle.
      */

      const auto maxMemoryAllocation =
          std::min(edgesComb,
                   globalBufferSize / uint64_t((sizeof(detail::DeviceTriplet) +
                                                sizeof(detail::SeedData)) *
                                               2));

      detail::DeviceTriplet* deviceCurvImpact =
          cl::sycl::malloc_device<detail::DeviceTriplet>(maxMemoryAllocation,
                                                         (*q));

      // Reserve memory in advance for seed indices and weight
      // Other way around would allocating it inside the loop
      // -> less memory usage, but more frequent allocation and deallocation
      detail::SeedData* deviceSeedArray =
          cl::sycl::malloc_device<detail::SeedData>(maxMemoryAllocation, *q);

      // Counting the seeds in the second kernel allows us to copy back the
      // right number of seeds, and no more.

      seeds.resize(M);
      uint32_t sumSeeds = 0;
      std::vector<uint32_t> deviceCountTriplets(edgesBottom, 0);

      // Do the triplet search and triplet filter for 2 sp fixed for middle
      // space points in the interval [firstMiddle, lastMiddle).

      uint32_t lastMiddle = 0;
      for (uint32_t firstMiddle = 0; firstMiddle < M;
           firstMiddle = lastMiddle) {
        // Determine the interval [firstMiddle, lastMiddle) right end based
        // on memory requirements.
        while (lastMiddle + 1 <= M && (sumBotTopCombined[lastMiddle + 1] -
                                           sumBotTopCombined[firstMiddle] <
                                       maxMemoryAllocation)) {
          ++lastMiddle;
        }

        const auto numTripletSearchThreads =
            sumBotTopCombined[lastMiddle] - sumBotTopCombined[firstMiddle];

        if (numTripletSearchThreads == 0)
          continue;

        sumSeeds = 0;
        deviceCountTriplets.resize(edgesBottom, 0);

        const auto numTripletFilterThreads =
            sumBotCompUptoMid[lastMiddle] - sumBotCompUptoMid[firstMiddle];

        const auto sumCombUptoFirstMiddle = sumBotTopCombined[firstMiddle];

        // Nd_range with maximum block size for triplet search and filter.
        // (global and local range is already given)
        cl::sycl::nd_range<1> tripletSearchNDRange =
            calculate1DimNDRange(numTripletSearchThreads, maxWorkGroupSize);

        cl::sycl::nd_range<1> tripletFilterNDRange =
            calculate1DimNDRange(numTripletFilterThreads, maxWorkGroupSize);

        sycl::buffer<uint32_t> countTripletsBuf(deviceCountTriplets.data(),
                                                edgesBottom);

        auto tripletKernel = q->submit([&](cl::sycl::handler& h) {
          h.depends_on({linB, linT});
          auto countTripletsAcc =
              sycl::ONEAPI::atomic_accessor<uint32_t, 1,
                                            sycl::ONEAPI::memory_order::relaxed,
                                            sycl::ONEAPI::memory_scope::device>(
                  countTripletsBuf, h);
          auto numTopDupletsAcc = numTopDupletsBuf.get_access<
              sycl::access::mode::read, sycl::access::target::global_buffer>(h);
          h.parallel_for<triplet_search_kernel>(
              tripletSearchNDRange, [=](cl::sycl::nd_item<1> item) {
                const uint32_t idx = item.get_global_linear_id();
                if (idx < numTripletSearchThreads) {
                  // Retrieve the index of the corresponding middle
                  // space point by binary search
                  auto L = firstMiddle;
                  auto R = lastMiddle;
                  auto mid = L;
                  while (L < R - 1) {
                    mid = (L + R) / 2;
                    // To be able to search in deviceSumComb, we need
                    // to use an offset (sumCombUptoFirstMiddle).
                    if (idx + sumCombUptoFirstMiddle < deviceSumComb[mid]) {
                      R = mid;
                    } else {
                      L = mid;
                    }
                  }
                  mid = L;

                  const auto numT = numTopDupletsAcc[mid];
                  const auto threadIdxForMiddleSP =
                      (idx - deviceSumComb[mid] + sumCombUptoFirstMiddle);

                  // NOTES ON THREAD MAPPING TO SPACE POINTS
                  /*
                    We need to map bottom and top SP indices to this
                    thread.

                    So we are mapping one bottom and one top SP to this thread
                    (we already have a middle SP) which gives us a tiplet.

                    This is done in the following way: We
                    calculated the number of possible triplet
                    combinations for this middle SP (let it be
                    num_comp_bot*num_comp_top). Let num_comp_bot = 2
                    and num_comp_top=3 in this example. So we have 2
                    compatible bottom and 3 compatible top SP for this
                    middle SP.

                    That gives us 6 threads altogether:
                               ===========================================
                    thread:    |  0   |  1   |  2   |  3   |  4   |  5   |
                    bottom id: | bot0 | bot0 | bot0 | bot1 | bot1 | bot1 |
                    top id:    | top0 | top1 | top2 | top0 | top1 | top2 |
                               ===========================================

                    If we divide 6 by the number of compatible top SP
                    for this middle SP, or deviceNumTopDuplets[mid]
                    which is 3 now, we get the id for the bottom SP.
                    Similarly, if we take modulo
                    deviceNumTopDuplets[mid], we get the id for the
                    top SP.

                    So if threadIdxForMiddleSP = 3, then ib = 1 and it = 0.

                    We can use these ids together with
                    deviceSumBot[mid] and deviceSumTop[mid] to be able
                    to index our other arrays.

                    These other arrays are deviceIndBot and deviceIndTop.

                    So to retrieve the bottom SP index for this thread, we'd
                    have to index the deviceIndBot array at
                      deviceSumBot[mid] + ib
                    which is the id for the bottom SP that we just calculated
                    (ib = 1 in the example).
                  */

                  const auto ib =
                      deviceSumBot[mid] + (threadIdxForMiddleSP / numT);
                  const auto it =
                      deviceSumTop[mid] + (threadIdxForMiddleSP % numT);

                  const auto linBotEq = deviceLinBot[ib];
                  const auto linTopEq = deviceLinTop[it];
                  const auto midSP = deviceMiddleSPs[mid];

                  const auto Vb = linBotEq.v;
                  const auto Ub = linBotEq.u;
                  const auto Erb = linBotEq.er;
                  const auto cotThetab = linBotEq.cotTheta;
                  const auto iDeltaRb = linBotEq.iDeltaR;

                  const auto Vt = linTopEq.v;
                  const auto Ut = linTopEq.u;
                  const auto Ert = linTopEq.er;
                  const auto cotThetat = linTopEq.cotTheta;
                  const auto iDeltaRt = linTopEq.iDeltaR;

                  const auto rM = midSP.r;
                  const auto varianceRM = midSP.varR;
                  const auto varianceZM = midSP.varZ;

                  auto iSinTheta2 = (1.f + cotThetab * cotThetab);
                  auto scatteringInRegion2 =
                      seedfinderConfig.maxScatteringAngle2 * iSinTheta2;
                  scatteringInRegion2 *= seedfinderConfig.sigmaScattering *
                                         seedfinderConfig.sigmaScattering;
                  auto error2 =
                      Ert + Erb +
                      2.f * (cotThetab * cotThetat * varianceRM + varianceZM) *
                          iDeltaRb * iDeltaRt;
                  auto deltaCotTheta = cotThetab - cotThetat;
                  auto deltaCotTheta2 = deltaCotTheta * deltaCotTheta;

                  deltaCotTheta = cl::sycl::abs(deltaCotTheta);
                  auto error = cl::sycl::sqrt(error2);
                  auto dCotThetaMinusError2 =
                      deltaCotTheta2 + error2 - 2.f * deltaCotTheta * error;
                  auto dU = Ut - Ub;

                  if ((!(deltaCotTheta2 - error2 > 0.f) ||
                       !(dCotThetaMinusError2 > scatteringInRegion2)) &&
                      !(dU == 0.f)) {
                    auto A = (Vt - Vb) / dU;
                    auto S2 = 1.f + A * A;
                    auto B = Vb - A * Ub;
                    auto B2 = B * B;

                    auto iHelixDiameter2 = B2 / S2;
                    auto pT2scatter =
                        4.f * iHelixDiameter2 * seedfinderConfig.pT2perRadius;
                    auto p2scatter = pT2scatter * iSinTheta2;
                    auto Im = cl::sycl::abs((A - B * rM) * rM);

                    if (!(S2 < B2 * seedfinderConfig.minHelixDiameter2) &&
                        !((deltaCotTheta2 - error2 > 0.f) &&
                          (dCotThetaMinusError2 >
                           p2scatter * seedfinderConfig.sigmaScattering *
                               seedfinderConfig.sigmaScattering)) &&
                        !(Im > seedfinderConfig.impactMax)) {
                      const auto top = deviceIndTop[it];
                      // this will be the t-th top space point for
                      // fixed middle and bottom SP
                      auto t = countTripletsAcc[ib].fetch_add(1);
                      /*
                        deviceSumComb[mid] - sumCombUptoFirstMiddle:
                        gives the memory location reserved for this
                        middle SP

                        (idx-deviceSumComb[mid]+sumCombUptoFirstMiddle:
                        this is the nth thread for this middle SP

                        (idx-deviceSumComb[mid]+sumCombUptoFirstMiddle)/numT:
                        this is the mth bottom SP for this middle SP

                        multiplying this by numT gives the memory
                        location for this middle and bottom SP

                        and by adding t to it, we will end up storing
                        compatible triplet candidates for this middle
                        and bottom SP right next to each other
                        starting from the given memory location
                      */
                      const auto tripletIdx = deviceSumComb[mid] -
                                              sumCombUptoFirstMiddle +
                                              (((idx - deviceSumComb[mid] +
                                                 sumCombUptoFirstMiddle) /
                                                numT) *
                                               numT) +
                                              t;

                      detail::DeviceTriplet T;
                      T.curvature = B / cl::sycl::sqrt(S2);
                      T.impact = Im;
                      T.topSPIndex = top;
                      deviceCurvImpact[tripletIdx] = T;
                    }
                  }
                }
              });
        });

        {
          sycl::buffer<uint32_t> countSeedsBuf(&sumSeeds, 1);
          q->submit([&](cl::sycl::handler& h) {
            h.depends_on(tripletKernel);
            auto countSeedsAcc = sycl::ONEAPI::atomic_accessor<
                uint32_t, 1, sycl::ONEAPI::memory_order::relaxed,
                sycl::ONEAPI::memory_scope::device>(countSeedsBuf, h);
            auto countTripletsAcc = countTripletsBuf.get_access<
                sycl::access::mode::read, sycl::access::target::global_buffer>(
                h);
            auto numTopDupletsAcc = numTopDupletsBuf.get_access<
                sycl::access::mode::read, sycl::access::target::global_buffer>(
                h);
            h.parallel_for<filter_2sp_fixed_kernel>(
                tripletFilterNDRange, [=](cl::sycl::nd_item<1> item) {
                  if (item.get_global_linear_id() < numTripletFilterThreads) {
                    const auto idx =
                        deviceSumBot[firstMiddle] + item.get_global_linear_id();
                    const auto mid = deviceMidIndPerBot[idx];
                    const auto bot = deviceIndBot[idx];

                    const auto tripletBegin =
                        deviceSumComb[mid] - sumCombUptoFirstMiddle +
                        (idx - deviceSumBot[mid]) * numTopDupletsAcc[mid];
                    const auto tripletEnd =
                        tripletBegin + countTripletsAcc[idx];

                    for (auto i1 = tripletBegin; i1 < tripletEnd; ++i1) {
                      const auto current = deviceCurvImpact[i1];
                      const auto top = current.topSPIndex;

                      const auto invHelixDiameter = current.curvature;
                      const auto lowerLimitCurv =
                          invHelixDiameter -
                          seedfinderConfig.deltaInvHelixDiameter;
                      const auto upperLimitCurv =
                          invHelixDiameter +
                          seedfinderConfig.deltaInvHelixDiameter;
                      const auto currentTop_r = deviceTopSPs[top].r;
                      auto weight = -(current.impact *
                                      seedfinderConfig.impactWeightFactor);

                      uint32_t compatCounter = 0;

                      // By default compatSeedLimit is 2 -> 2 is
                      // currently hard coded, because variable length
                      // arrays are not supported in SYCL kernels.
                      float compatibleSeedR[2];
                      for (auto i2 = tripletBegin;
                           i2 < tripletEnd &&
                           compatCounter < seedfinderConfig.compatSeedLimit;
                           ++i2) {
                        const auto other = deviceCurvImpact[i2];

                        const auto otherCurv = other.curvature;
                        const auto otherTop_r =
                            deviceTopSPs[other.topSPIndex].r;
                        const float deltaR =
                            cl::sycl::abs(currentTop_r - otherTop_r);
                        if (deltaR >= seedfinderConfig.filterDeltaRMin &&
                            otherCurv >= lowerLimitCurv &&
                            otherCurv <= upperLimitCurv) {
                          uint32_t c = 0;
                          for (;
                               c < compatCounter &&
                               cl::sycl::abs(compatibleSeedR[c] - otherTop_r) >=
                                   seedfinderConfig.filterDeltaRMin;
                               ++c) {
                          }
                          if (c == compatCounter) {
                            compatibleSeedR[c] = otherTop_r;
                            ++compatCounter;
                          }
                        }
                      }

                      weight +=
                          compatCounter * seedfinderConfig.compatSeedWeight;

                      const auto bottomSP = deviceBottomSPs[bot];
                      const auto middleSP = deviceMiddleSPs[mid];
                      const auto topSP = deviceTopSPs[top];

                      weight +=
                          deviceCuts.seedWeight(bottomSP, middleSP, topSP);

                      if (deviceCuts.singleSeedCut(weight, bottomSP, middleSP,
                                                   topSP)) {
                        const auto i = countSeedsAcc[0].fetch_add(1);
                        detail::SeedData D;
                        D.bottom = bot;
                        D.top = top;
                        D.middle = mid;
                        D.weight = weight;
                        deviceSeedArray[i] = D;
                      }
                    }
                  }
                });
          });
        }  // sync, countSeedsBuf gets destroyed and its value is copied back to
           // sumSeeds

        if (sumSeeds != 0) {
          std::vector<detail::SeedData> hostSeedArray(sumSeeds);
          auto e0 = q->memcpy(&hostSeedArray[0], deviceSeedArray,
                              sumSeeds * sizeof(detail::SeedData));
          e0.wait();

          for (uint32_t t = 0; t < sumSeeds; ++t) {
            auto m = hostSeedArray[t].middle;
            seeds[m].push_back(hostSeedArray[t]);
          }
        }
      }

      //************************************************//
      // ************ TRIPLET SEARCH - END ************ //
      //************************************************//

      // NOTES ON MEMORY MANAGEMENT
      /*
        Note that memory management is very naive, but this should only be a
        temporary solution. A better idea is to use a separate memory
        manager, see the CUDA (2) implementation.
      */

      cl::sycl::free(deviceLinBot, *q);
      cl::sycl::free(deviceLinTop, *q);

      cl::sycl::free(deviceIndBot, *q);
      cl::sycl::free(deviceIndTop, *q);
      cl::sycl::free(deviceMidIndPerBot, *q);
      cl::sycl::free(deviceMidIndPerTop, *q);
      cl::sycl::free(deviceSumBot, *q);
      cl::sycl::free(deviceSumTop, *q);
      cl::sycl::free(deviceSumComb, *q);

      cl::sycl::free(deviceCurvImpact, *q);
      cl::sycl::free(deviceSeedArray, *q);
    }

    cl::sycl::free(deviceTmpIndBot, *q);
    cl::sycl::free(deviceTmpIndTop, *q);
    cl::sycl::free(deviceBottomSPs, *q);
    cl::sycl::free(deviceMiddleSPs, *q);
    cl::sycl::free(deviceTopSPs, *q);
  } catch (cl::sycl::exception const& e) {
    ACTS_LOCAL_LOGGER(
        Acts::getDefaultLogger("SyclSeeding", Acts::Logging::INFO));
    ACTS_FATAL("Caught synchronous SYCL exception:\n" << e.what())
    exit(0);
  }
};
}  // namespace Acts::Sycl

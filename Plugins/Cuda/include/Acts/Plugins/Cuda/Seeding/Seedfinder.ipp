// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <cmath>
#include <numeric>
#include <type_traits>

namespace Acts {

template <typename external_spacepoint_t>
Seedfinder<external_spacepoint_t, Acts::Cuda>::Seedfinder(
    Acts::SeedfinderConfig<external_spacepoint_t> config)
    : m_config(std::move(config)) {
  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  m_config.highland = 13.6 * std::sqrt(m_config.radLengthPerSeed) *
                      (1 + 0.038 * std::log(m_config.radLengthPerSeed));
  float maxScatteringAngle = m_config.highland / m_config.minPt;
  m_config.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;

  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
  // millimeter
  // TODO: change using ACTS units
  m_config.pTPerHelixRadius = 300. * m_config.bFieldInZ;
  m_config.minHelixDiameter2 =
      std::pow(m_config.minPt * 2 / m_config.pTPerHelixRadius, 2);
  m_config.pT2perRadius =
      std::pow(m_config.highland / m_config.pTPerHelixRadius, 2);
}

// CUDA seed finding
template <typename external_spacepoint_t>
template <typename sp_range_t>
std::vector<Seed<external_spacepoint_t>>
Seedfinder<external_spacepoint_t, Acts::Cuda>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const {
  std::vector<Seed<external_spacepoint_t>> outputVec;

  // Get SeedfinderConfig values
  CudaScalar<float> deltaRMin_cuda(&m_config.deltaRMin);
  CudaScalar<float> deltaRMax_cuda(&m_config.deltaRMax);
  CudaScalar<float> cotThetaMax_cuda(&m_config.cotThetaMax);
  CudaScalar<float> collisionRegionMin_cuda(&m_config.collisionRegionMin);
  CudaScalar<float> collisionRegionMax_cuda(&m_config.collisionRegionMax);
  CudaScalar<float> maxScatteringAngle2_cuda(&m_config.maxScatteringAngle2);
  CudaScalar<float> sigmaScattering_cuda(&m_config.sigmaScattering);
  CudaScalar<float> minHelixDiameter2_cuda(&m_config.minHelixDiameter2);
  CudaScalar<float> pT2perRadius_cuda(&m_config.pT2perRadius);
  CudaScalar<float> impactMax_cuda(&m_config.impactMax);
  auto seedFilterConfig = m_config.seedFilter->m_cfg;  
  CudaScalar<float> deltaInvHelixDiameter_cuda(&seedFilterConfig.deltaInvHelixDiameter);
  CudaScalar<float> impactWeightFactor_cuda(&seedFilterConfig.impactWeightFactor);
  
  //---------------------------------
  // Algorithm 0. Matrix Flattening
  //---------------------------------

  std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*>
      middleSPvec;
  std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*>
      bottomSPvec;
  std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*> topSPvec;

  // Get the size of spacepoints
  int nSpM(0);
  int nSpB(0);
  int nSpT(0);

  for (auto sp : middleSPs) {
    nSpM++;
    middleSPvec.push_back(sp);
  }
  for (auto sp : bottomSPs) {
    nSpB++;
    bottomSPvec.push_back(sp);
  }
  for (auto sp : topSPs) {
    nSpT++;
    topSPvec.push_back(sp);
  }

  CudaScalar<int> nSpM_cuda(&nSpM);
  CudaScalar<int> nSpB_cuda(&nSpB);
  CudaScalar<int> nSpT_cuda(&nSpT);

  if (nSpM == 0 || nSpB == 0 || nSpT == 0)
    return outputVec;

  // Matrix flattening
  CpuMatrix<float> spMmat_cpu(nSpM, 6);  // x y z r varR varZ
  CpuMatrix<float> spBmat_cpu(nSpB, 6);
  CpuMatrix<float> spTmat_cpu(nSpT, 6);

  auto fillMatrix = [](auto& mat, auto& id, auto sp) {
    mat.set(id, 0, sp->x());
    mat.set(id, 1, sp->y());
    mat.set(id, 2, sp->z());
    mat.set(id, 3, sp->radius());
    mat.set(id, 4, sp->varianceR());
    mat.set(id, 5, sp->varianceZ());
    id++;
  };

  int mIdx(0);
  int bIdx(0);
  int tIdx(0);

  for (auto sp : middleSPs) {
    fillMatrix(spMmat_cpu, mIdx, sp);
  }

  for (auto sp : bottomSPs) {
    fillMatrix(spBmat_cpu, bIdx, sp);
  }

  for (auto sp : topSPs) {
    fillMatrix(spTmat_cpu, tIdx, sp);
  }

  CudaMatrix<float> spMmat_cuda(nSpM, 6, &spMmat_cpu);
  CudaMatrix<float> spBmat_cuda(nSpB, 6, &spBmat_cpu);
  CudaMatrix<float> spTmat_cuda(nSpT, 6, &spTmat_cpu);

  //------------------------------------
  //  Algorithm 1. Doublet Search (DS)
  //------------------------------------

  CudaScalar<int> nSpMcomp_cuda(new int(0));
  CudaScalar<int> nSpBcompPerSpM_Max_cuda(new int(0));
  CudaScalar<int> nSpTcompPerSpM_Max_cuda(new int(0));
  CudaVector<int> nSpBcompPerSpM_cuda(nSpM);
  nSpBcompPerSpM_cuda.zeros();
  CudaVector<int> nSpTcompPerSpM_cuda(nSpM);
  nSpTcompPerSpM_cuda.zeros();
  CudaVector<int> McompIndex_cuda(nSpM);
  CudaMatrix<int> BcompIndex_cuda(nSpB, nSpM);
  CudaMatrix<int> TcompIndex_cuda(nSpT, nSpM);
  CudaMatrix<int> tmpBcompIndex_cuda(nSpB, nSpM);
  CudaMatrix<int> tmpTcompIndex_cuda(nSpT, nSpM);

  dim3 DS_BlockSize = m_config.maxBlockSize;
  dim3 DS_GridSize(nSpM, 1, 1);

  searchDoublet(DS_GridSize, DS_BlockSize, nSpM_cuda.get(), spMmat_cuda.get(),
                nSpB_cuda.get(), spBmat_cuda.get(), nSpT_cuda.get(),
                spTmat_cuda.get(), deltaRMin_cuda.get(), deltaRMax_cuda.get(),
                cotThetaMax_cuda.get(), collisionRegionMin_cuda.get(),
                collisionRegionMax_cuda.get(), nSpMcomp_cuda.get(),
                nSpBcompPerSpM_Max_cuda.get(), nSpTcompPerSpM_Max_cuda.get(),
                nSpBcompPerSpM_cuda.get(), nSpTcompPerSpM_cuda.get(),
                McompIndex_cuda.get(), BcompIndex_cuda.get(),
                tmpBcompIndex_cuda.get(), TcompIndex_cuda.get(),
                tmpTcompIndex_cuda.get());

  CpuScalar<int> nSpMcomp_cpu(&nSpMcomp_cuda);
  CpuScalar<int> nSpBcompPerSpM_Max_cpu(&nSpBcompPerSpM_Max_cuda);
  CpuScalar<int> nSpTcompPerSpM_Max_cpu(&nSpTcompPerSpM_Max_cuda);
  CpuVector<int> nSpBcompPerSpM_cpu(nSpM, &nSpBcompPerSpM_cuda);
  CpuVector<int> nSpTcompPerSpM_cpu(nSpM, &nSpTcompPerSpM_cuda);
  CpuVector<int> McompIndex_cpu(nSpM, &McompIndex_cuda);
  CpuMatrix<int> BcompIndex_cpu(nSpB, nSpM, &BcompIndex_cuda);
  CpuMatrix<int> TcompIndex_cpu(nSpT, nSpM, &TcompIndex_cuda);

  //--------------------------------------
  //  Algorithm 2. Transform coordinate
  //--------------------------------------

  CudaMatrix<float> spMcompMat_cuda(*nSpMcomp_cpu.get(), 6);

  CudaMatrix<float> spBcompMatPerSpM_cuda(*nSpBcompPerSpM_Max_cpu.get(),
                                          (*nSpMcomp_cpu.get()) * 6);
  CudaMatrix<float> spTcompMatPerSpM_cuda(*nSpTcompPerSpM_Max_cpu.get(),
                                          (*nSpMcomp_cpu.get()) * 6);
  CudaMatrix<float> circBcompMatPerSpM_cuda(*nSpBcompPerSpM_Max_cpu.get(),
                                            (*nSpMcomp_cpu.get()) * 6);
  CudaMatrix<float> circTcompMatPerSpM_cuda(*nSpTcompPerSpM_Max_cpu.get(),
                                            (*nSpMcomp_cpu.get()) * 6);

  dim3 TC_GridSize(*nSpMcomp_cpu.get(), 1, 1);
  dim3 TC_BlockSize = m_config.maxBlockSize;

  transformCoordinate(
      TC_GridSize, TC_BlockSize, nSpM_cuda.get(), spMmat_cuda.get(),
      McompIndex_cuda.get(), nSpB_cuda.get(), spBmat_cuda.get(),
      nSpBcompPerSpM_Max_cuda.get(), BcompIndex_cuda.get(), nSpT_cuda.get(),
      spTmat_cuda.get(), nSpTcompPerSpM_Max_cuda.get(), TcompIndex_cuda.get(),
      spMcompMat_cuda.get(), spBcompMatPerSpM_cuda.get(),
      circBcompMatPerSpM_cuda.get(), spTcompMatPerSpM_cuda.get(),
      circTcompMatPerSpM_cuda.get());

  //------------------------------------------------------
  //  Algorithm 3. Triplet Search (TS) & Seed filtering
  //------------------------------------------------------

  // retreive middle-bottom doublet circ information
  CpuMatrix<float> circBcompMatPerSpM_cpu(*nSpBcompPerSpM_Max_cpu.get(),
                                          (*nSpMcomp_cpu.get()) * 6,
                                          &circBcompMatPerSpM_cuda);

  std::vector<int> offsetVec(100);
  std::iota(std::begin(offsetVec), std::end(offsetVec),
            0);  // Fill with 0, 1, ..., 99.
  for (auto& el : offsetVec)
    el = el * m_config.maxBlockSize;
  CudaVector<int> offsetVec_cuda(offsetVec.size(), &offsetVec[0]);

  const int nTrplPerSpMLimit =
      m_config.nTrplPerSpBLimit * (*nSpBcompPerSpM_Max_cpu.get());
  CudaScalar<int> nTrplPerSpMLimit_cuda(&nTrplPerSpMLimit);

  CudaVector<int> nTrplPerSpM_cuda(*nSpMcomp_cpu.get());
  nTrplPerSpM_cuda.zeros();
  CudaMatrix<int> TtrplIndex_cuda(nTrplPerSpMLimit, *nSpMcomp_cpu.get());
  CudaMatrix<int> BtrplIndex_cuda(nTrplPerSpMLimit, *nSpMcomp_cpu.get());
  CudaMatrix<float> curvatures_cuda(nTrplPerSpMLimit, *nSpMcomp_cpu.get());
  CudaMatrix<float> impactparameters_cuda(nTrplPerSpMLimit,
                                          *nSpMcomp_cpu.get());

  CpuVector<int> nTrplPerSpM_cpu(*nSpMcomp_cpu.get(), true);
  nTrplPerSpM_cpu.zeros();
  CpuMatrix<int> TtrplIndex_cpu(nTrplPerSpMLimit, *nSpMcomp_cpu.get(), true);
  CpuMatrix<int> BtrplIndex_cpu(nTrplPerSpMLimit, *nSpMcomp_cpu.get(), true);
  CpuMatrix<float> curvatures_cpu(nTrplPerSpMLimit, *nSpMcomp_cpu.get(), true);
  CpuMatrix<float> impactparameters_cpu(nTrplPerSpMLimit, *nSpMcomp_cpu.get(),
                                        true);

  cudaStream_t cuStream;
  cudaStreamCreate(&cuStream);

  for (int i_m = 0; i_m <= *nSpMcomp_cpu.get(); i_m++) {
    cudaStreamSynchronize(cuStream);

    // Search Triplet
    if (i_m < *nSpMcomp_cpu.get()) {
      int mIndex = *McompIndex_cpu.get(i_m);
      int nSpBcompPerSpM = *nSpBcompPerSpM_cpu.get(mIndex);
      int nSpTcompPerSpM = *nSpTcompPerSpM_cpu.get(mIndex);

      dim3 TS_GridSize(nSpBcompPerSpM, 1, 1);
      dim3 TS_BlockSize;

      int i_ts(0);
      while (offsetVec[i_ts] < nSpTcompPerSpM) {
        TS_BlockSize =
            dim3(fmin(m_config.maxBlockSize, nSpTcompPerSpM - offsetVec[i_ts]),
                 1, 1);
        searchTriplet(TS_GridSize, TS_BlockSize, offsetVec_cuda.get(i_ts),
                      nSpMcomp_cuda.get(), spMcompMat_cuda.get(i_m, 0),
                      nSpBcompPerSpM_Max_cuda.get(),
                      circBcompMatPerSpM_cuda.get(0, 6 * i_m),
                      nSpTcompPerSpM_Max_cuda.get(),
		      spTcompMatPerSpM_cuda.get(0, 6 * i_m),
                      circTcompMatPerSpM_cuda.get(0, 6 * i_m),
                      // Seed finder config
                      maxScatteringAngle2_cuda.get(),
                      sigmaScattering_cuda.get(), minHelixDiameter2_cuda.get(),
                      pT2perRadius_cuda.get(), impactMax_cuda.get(),
                      nTrplPerSpMLimit_cuda.get(),
		      deltaInvHelixDiameter_cuda.get(),
		      impactWeightFactor_cuda.get(),		      
                      // output
                      nTrplPerSpM_cuda.get(i_m), TtrplIndex_cuda.get(0, i_m),
                      BtrplIndex_cuda.get(0, i_m), curvatures_cuda.get(0, i_m),
                      impactparameters_cuda.get(0, i_m), &cuStream);
        i_ts++;
      }

      nTrplPerSpM_cpu.copyD2H(nTrplPerSpM_cuda.get(i_m), 1, i_m, &cuStream);

      TtrplIndex_cpu.copyD2H(TtrplIndex_cuda.get(0, i_m), nTrplPerSpMLimit,
                             nTrplPerSpMLimit * i_m, &cuStream);

      BtrplIndex_cpu.copyD2H(BtrplIndex_cuda.get(0, i_m), nTrplPerSpMLimit,
                             nTrplPerSpMLimit * i_m, &cuStream);

      curvatures_cpu.copyD2H(curvatures_cuda.get(0, i_m), nTrplPerSpMLimit,
                             nTrplPerSpMLimit * i_m, &cuStream);

      impactparameters_cpu.copyD2H(impactparameters_cuda.get(0, i_m),
                                   nTrplPerSpMLimit, nTrplPerSpMLimit * i_m,
                                   &cuStream);
    }

    std::vector<std::pair<
        float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
        seedsPerSpM;

    // Seed filtering
    if (i_m > 0) {
      int i_prev = i_m - 1;
      auto mIndex = *McompIndex_cpu.get(i_prev);
      auto bIndexVec = BcompIndex_cpu.get(0, i_prev);
      auto tIndexVec = TcompIndex_cpu.get(0, i_prev);

      // Bottom index (key) is sorted in the map automatrically
      std::map<int, std::vector<std::tuple<int, float, float>>> trplMap;

      for (int i_trpl = 0; i_trpl < *nTrplPerSpM_cpu.get(i_prev); i_trpl++) {
        int BtrplIndex = *(BtrplIndex_cpu.get(i_trpl, i_prev));
        int TriplIndex = *(TtrplIndex_cpu.get(i_trpl, i_prev));
        int bIndex = bIndexVec[BtrplIndex];
        int tIndex = tIndexVec[TriplIndex];

        float curv = *(curvatures_cpu.get(i_trpl, i_prev));
        float impact = *(impactparameters_cpu.get(i_trpl, i_prev));

        trplMap[bIndex].push_back(std::make_tuple(tIndex, curv, impact));
      }

      int nTrplPerSpBMax(0);
      // Sort top index
      for (auto& trplPerSpB : trplMap) {
        sort(trplPerSpB.second.begin(), trplPerSpB.second.end());
        int nTrplPerSpB = trplPerSpB.second.size();
        nTrplPerSpBMax = std::max(nTrplPerSpB, nTrplPerSpBMax);
      }

      std::vector<const InternalSpacePoint<external_spacepoint_t>*> tVec;
      std::vector<float> curvatures;
      std::vector<float> impactParameters;
      tVec.reserve(nTrplPerSpBMax);
      curvatures.reserve(nTrplPerSpBMax);
      impactParameters.reserve(nTrplPerSpBMax);

      for (auto& trplPerSpB : trplMap) {
        tVec.clear();
        curvatures.clear();
        impactParameters.clear();

        int bIndex = trplPerSpB.first;
        auto triplets = trplPerSpB.second;

        for (auto it = std::make_move_iterator(triplets.begin()),
                  end = std::make_move_iterator(triplets.end());
             it != end; ++it) {
          int tIndex = std::get<0>(*it);
          tVec.push_back(std::move(topSPvec[tIndex]));
          curvatures.push_back(std::move(std::get<1>(*it)));
          impactParameters.push_back(std::move(std::get<2>(*it)));
        }

        float Zob = *(circBcompMatPerSpM_cpu.get(bIndex, (i_prev)*6));

        std::vector<std::pair<
            float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
            sameTrackSeeds;
        sameTrackSeeds = std::move(m_config.seedFilter->filterSeeds_2SpFixed(
            *bottomSPvec[bIndex], *middleSPvec[mIndex], tVec, curvatures,
            impactParameters, Zob));
        seedsPerSpM.insert(seedsPerSpM.end(),
                           std::make_move_iterator(sameTrackSeeds.begin()),
                           std::make_move_iterator(sameTrackSeeds.end()));
      }
      m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, outputVec);
    }
  }
  return outputVec;
}
}  // namespace Acts

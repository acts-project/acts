// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Cuda.hpp"
#include "Acts/Plugins/Cuda/Seeding/Kernels.cuh"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <type_traits>

namespace Acts {
namespace Cuda {

template <typename external_spacepoint_t>
Seedfinder<external_spacepoint_t>::Seedfinder(
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
Seedfinder<external_spacepoint_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const {
  std::vector<Seed<external_spacepoint_t>> outputVec;

  // Get SeedfinderConfig values
  const auto seedFilterConfig = m_config.seedFilter->getSeedFilterConfig();

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

  if (nSpM == 0 || nSpB == 0 || nSpT == 0)
    return outputVec;

  // Matrix flattening
  HostMatrix<float> spMmat_cpu(nSpM, 6);  // x y z r varR varZ
  HostMatrix<float> spBmat_cpu(nSpB, 6);
  HostMatrix<float> spTmat_cpu(nSpT, 6);

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
  for (auto sp : middleSPs) {
    fillMatrix(spMmat_cpu, mIdx, sp);
  }
  int bIdx(0);
  for (auto sp : bottomSPs) {
    fillMatrix(spBmat_cpu, bIdx, sp);
  }
  int tIdx(0);
  for (auto sp : topSPs) {
    fillMatrix(spTmat_cpu, tIdx, sp);
  }

  DeviceMatrix<float> spMmat_cuda(nSpM, 6);
  copyToDevice(spMmat_cuda, spMmat_cpu);
  DeviceMatrix<float> spBmat_cuda(nSpB, 6);
  copyToDevice(spBmat_cuda, spBmat_cpu);
  DeviceMatrix<float> spTmat_cuda(nSpT, 6);
  copyToDevice(spTmat_cuda, spTmat_cpu);
  //------------------------------------
  //  Algorithm 1. Doublet Search (DS)
  //------------------------------------

  ResultScalar<int> nSpMcomp_cuda;
  ResultScalar<int> nSpBcompPerSpMMax_cuda;
  ResultScalar<int> nSpTcompPerSpMMax_cuda;
  DeviceVector<int> nSpBcompPerSpM_cuda(nSpM);
  nSpBcompPerSpM_cuda.zeros();
  DeviceVector<int> nSpTcompPerSpM_cuda(nSpM);
  nSpTcompPerSpM_cuda.zeros();
  DeviceVector<int> McompIndex_cuda(nSpM);
  DeviceMatrix<int> BcompIndex_cuda(nSpB, nSpM);
  DeviceMatrix<int> TcompIndex_cuda(nSpT, nSpM);
  DeviceMatrix<int> tmpBcompIndex_cuda(nSpB, nSpM);
  DeviceMatrix<int> tmpTcompIndex_cuda(nSpT, nSpM);

  dim3 DS_BlockSize = m_config.maxBlockSize;
  dim3 DS_GridSize(nSpM, 1, 1);

  searchDoublet(DS_GridSize, DS_BlockSize,
                nSpM, spMmat_cuda.getPtr(),
                nSpB, spBmat_cuda.getPtr(),
                nSpT, spTmat_cuda.getPtr(),
                m_config.deltaRMin, m_config.deltaRMax,
                m_config.cotThetaMax, m_config.collisionRegionMin,
                m_config.collisionRegionMax,
                nSpMcomp_cuda.getPtr(),
                nSpBcompPerSpMMax_cuda.getPtr(), nSpTcompPerSpMMax_cuda.getPtr(),
                nSpBcompPerSpM_cuda.getPtr(), nSpTcompPerSpM_cuda.getPtr(),
                McompIndex_cuda.getPtr(), BcompIndex_cuda.getPtr(),
                tmpBcompIndex_cuda.getPtr(), TcompIndex_cuda.getPtr(),
                tmpTcompIndex_cuda.getPtr());

  const int nSpMcomp_cpu = nSpMcomp_cuda;
  const int nSpBcompPerSpMMax_cpu = nSpBcompPerSpMMax_cuda;
  const int nSpTcompPerSpMMax_cpu = nSpTcompPerSpMMax_cuda;
  HostVector<int> nSpBcompPerSpM_cpu(nSpM);
  copyToHost(nSpBcompPerSpM_cpu, nSpBcompPerSpM_cuda);
  HostVector<int> nSpTcompPerSpM_cpu(nSpM);
  copyToHost(nSpTcompPerSpM_cpu, nSpTcompPerSpM_cuda);
  HostVector<int> McompIndex_cpu(nSpM);
  copyToHost(McompIndex_cpu, McompIndex_cuda);

  //--------------------------------------
  //  Algorithm 2. Transform coordinate
  //--------------------------------------

  DeviceMatrix<float> spMcompMat_cuda(nSpMcomp_cpu, 6);
  DeviceMatrix<float> spBcompMatPerSpM_cuda(nSpBcompPerSpMMax_cpu,
                                            nSpMcomp_cpu * 6);
  DeviceMatrix<float> spTcompMatPerSpM_cuda(nSpTcompPerSpMMax_cpu,
                                            nSpMcomp_cpu * 6);
  DeviceMatrix<float> circBcompMatPerSpM_cuda(nSpBcompPerSpMMax_cpu,
                                              nSpMcomp_cpu * 6);
  DeviceMatrix<float> circTcompMatPerSpM_cuda(nSpTcompPerSpMMax_cpu,
                                              nSpMcomp_cpu * 6);

  dim3 TC_GridSize(nSpMcomp_cpu, 1, 1);
  dim3 TC_BlockSize = m_config.maxBlockSize;

  transformCoordinate(
      TC_GridSize, TC_BlockSize, nSpM, spMmat_cuda.getPtr(),
      McompIndex_cuda.getPtr(), nSpB, spBmat_cuda.getPtr(),
      nSpBcompPerSpMMax_cpu, BcompIndex_cuda.getPtr(), nSpT,
      spTmat_cuda.getPtr(), nSpTcompPerSpMMax_cpu, TcompIndex_cuda.getPtr(),
      spMcompMat_cuda.getPtr(), spBcompMatPerSpM_cuda.getPtr(),
      circBcompMatPerSpM_cuda.getPtr(), spTcompMatPerSpM_cuda.getPtr(),
      circTcompMatPerSpM_cuda.getPtr());

  //------------------------------------------------------
  //  Algorithm 3. Triplet Search (TS) & Seed filtering
  //------------------------------------------------------

  const int nTrplPerSpMLimit =
      m_config.nAvgTrplPerSpBLimit * nSpBcompPerSpMMax_cpu;

  DeviceVector<int> nTrplPerSpM_cuda(nSpMcomp_cpu);
  nTrplPerSpM_cuda.zeros();
  DeviceMatrix<Triplet> TripletsPerSpM_cuda(nTrplPerSpMLimit,
                                            nSpMcomp_cpu);
  HostVector<int> nTrplPerSpM_cpu(nSpMcomp_cpu);
  nTrplPerSpM_cpu.zeros();
  HostMatrix<Triplet> TripletsPerSpM_cpu(nTrplPerSpMLimit, nSpMcomp_cpu);
  cudaStream_t cuStream;
  cudaStreamCreate(&cuStream);

  for (int i_m = 0; i_m <= nSpMcomp_cpu; i_m++) {
    cudaStreamSynchronize(cuStream);

    // Search Triplet
    if (i_m < nSpMcomp_cpu) {
      const int mIndex = McompIndex_cpu.get(i_m);
      const int nSpBcompPerSpM = nSpBcompPerSpM_cpu.get(mIndex);
      const int nSpTcompPerSpM = nSpTcompPerSpM_cpu.get(mIndex);

      dim3 TS_GridSize(nSpBcompPerSpM, 1, 1);
      dim3 TS_BlockSize =
          dim3(fmin(m_config.maxBlockSize, nSpTcompPerSpM), 1, 1);

      searchTriplet(
          TS_GridSize, TS_BlockSize,
          nSpTcompPerSpM_cuda.getPtr(mIndex), nSpMcomp_cpu,
          spMcompMat_cuda.getPtr(i_m, 0), nSpBcompPerSpMMax_cpu,
          BcompIndex_cuda.getPtr(0, i_m), circBcompMatPerSpM_cuda.getPtr(0, 6 * i_m),
          nSpTcompPerSpMMax_cpu, TcompIndex_cuda.getPtr(0, i_m),
          spTcompMatPerSpM_cuda.getPtr(0, 6 * i_m),
          circTcompMatPerSpM_cuda.getPtr(0, 6 * i_m),
          // Seed finder config
          m_config.maxScatteringAngle2, m_config.sigmaScattering,
          m_config.minHelixDiameter2, m_config.pT2perRadius,
          m_config.impactMax, nTrplPerSpMLimit,
          m_config.nTrplPerSpBLimit,
          seedFilterConfig.deltaInvHelixDiameter, seedFilterConfig.impactWeightFactor,
          seedFilterConfig.deltaRMin, seedFilterConfig.compatSeedWeight,
          seedFilterConfig.compatSeedLimit,
          // output
          nTrplPerSpM_cuda.getPtr(i_m), TripletsPerSpM_cuda.getPtr(0, i_m),
          cuStream);
      nTrplPerSpM_cpu.copyFrom(nTrplPerSpM_cuda.getPtr(i_m), 1, i_m, cuStream);

      TripletsPerSpM_cpu.copyFrom(TripletsPerSpM_cuda.getPtr(0, i_m),
                                  nTrplPerSpMLimit, nTrplPerSpMLimit * i_m,
                                  cuStream);
    }

    if (i_m > 0) {
      const auto m_experimentCuts = m_config.seedFilter->getExperimentCuts();
      std::vector<std::pair<
          float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
          seedsPerSpM;

      for (int i = 0; i < nTrplPerSpM_cpu.get(i_m - 1); i++) {
        auto& triplet = TripletsPerSpM_cpu.get(i, i_m - 1);
        int mIndex = McompIndex_cpu.get(i_m - 1);
        int bIndex = triplet.bIndex;
        int tIndex = triplet.tIndex;

        auto& bottomSP = *bottomSPvec[bIndex];
        auto& middleSP = *middleSPvec[mIndex];
        auto& topSP = *topSPvec[tIndex];
        if (m_experimentCuts != nullptr) {
          // add detector specific considerations on the seed weight
          triplet.weight +=
              m_experimentCuts->seedWeight(bottomSP, middleSP, topSP);
          // discard seeds according to detector specific cuts (e.g.: weight)
          if (!m_experimentCuts->singleSeedCut(triplet.weight, bottomSP,
                                               middleSP, topSP)) {
            continue;
          }
        }

        float Zob = 0;  // It is not used in the seed filter but needs to be
                        // fixed anyway...

        seedsPerSpM.push_back(std::make_pair(
            triplet.weight,
            std::make_unique<const InternalSeed<external_spacepoint_t>>(
                bottomSP, middleSP, topSP, Zob)));
      }

      m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, outputVec);
    }
  }
  return outputVec;
}

} // namespace Cuda
} // namespace Acts

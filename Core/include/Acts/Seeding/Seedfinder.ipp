// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include <numeric>
#include <type_traits>
#include <algorithm>
#include <chrono>
#include <numeric>

namespace Acts {

  template <typename external_spacepoint_t, typename platform_t>
  Seedfinder<external_spacepoint_t, platform_t>::Seedfinder(
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

  t_metric = std::make_tuple(0,0,0,0);
  
  }

  
  template< typename external_spacepoint_t, typename platform_t>
  template< typename T, typename sp_range_t>
  typename std::enable_if< std::is_same<T, Acts::CPU>::value, std::vector<Seed<external_spacepoint_t> > >::type
  Seedfinder<external_spacepoint_t, platform_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const{
  std::vector<Seed<external_spacepoint_t>> outputVec;
  auto start_wall = std::chrono::system_clock::now();

  //int BcompMax=0;
  //int TcompMax=0;
  
  for (auto spM : middleSPs) {    

    auto start_DS = std::chrono::system_clock::now();
    
    float rM = spM->radius();
    float zM = spM->z();
    float varianceRM = spM->varianceR();
    float varianceZM = spM->varianceZ();
    
    // Doublet search    
    auto compatBottomSP =
      SeedfinderCpuFunctions<external_spacepoint_t,
			     sp_range_t>::searchDoublet(true, bottomSPs, *spM, m_config);
    
    // no bottom SP found -> try next spM
    if (compatBottomSP.empty()) {
      continue;
    }

    auto compatTopSP =
      SeedfinderCpuFunctions<external_spacepoint_t,
			     sp_range_t>::searchDoublet(false, topSPs, *spM, m_config);

    // no top SP found -> try next spM
    if (compatTopSP.empty()) {
      continue;
    }
    
    auto end_DS = std::chrono::system_clock::now();
    std::chrono::duration<double> elapse_DS = end_DS-start_DS;
    std::get<0>(t_metric) += elapse_DS.count();
    
    // contains parameters required to calculate circle with linear equation

    auto start_TC = std::chrono::system_clock::now();
    
    // ...for bottom-middle
    std::vector<LinCircle> linCircleBottom;
    // ...for middle-top
    std::vector<LinCircle> linCircleTop;
    
    SeedfinderCpuFunctions<external_spacepoint_t,sp_range_t>::transformCoordinates(compatBottomSP, *spM, true, linCircleBottom);
    SeedfinderCpuFunctions<external_spacepoint_t,sp_range_t>::transformCoordinates(compatTopSP, *spM, false, linCircleTop);

    auto end_TC = std::chrono::system_clock::now();
    std::chrono::duration<double> elapse_TC = end_TC-start_TC;
    std::get<1>(t_metric) += elapse_TC.count();

    auto start_TS = std::chrono::system_clock::now();
    
    auto seedsPerSpM = SeedfinderCpuFunctions<external_spacepoint_t,sp_range_t>::searchTriplet(*spM, compatBottomSP, compatTopSP, linCircleBottom, linCircleTop, m_config);
    m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, outputVec);

    auto end_TS = std::chrono::system_clock::now();
    std::chrono::duration<double> elapse_TS = end_TS-start_TS;
    std::get<2>(t_metric) += (elapse_TS).count();
  }  

  auto end_wall = std::chrono::system_clock::now();
  std::chrono::duration<double> elapse_wall = end_wall-start_wall;
  std::get<3>(t_metric) += elapse_wall.count();
  
  return outputVec;
  }

#ifdef  ACTS_HAS_CUDA
  
  // CUDA seed finding
  template< typename external_spacepoint_t, typename platform_t>
  template< typename T, typename sp_range_t>
  typename std::enable_if< std::is_same<T, Acts::CUDA>::value, std::vector<Seed<external_spacepoint_t> > >::type
  Seedfinder<external_spacepoint_t, platform_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const{
  std::vector<Seed<external_spacepoint_t>> outputVec;

  // Get SeedfinderConfig values
  CudaScalar<bool>  true_cuda(new bool(true));  
  CudaScalar<bool>  false_cuda(new bool(false));    
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
  
  /*----------------------------------
     Algorithm 0. Matrix Flattening 
  ----------------------------------*/
  
  auto start_wall = std::chrono::system_clock::now();
  auto start_DS = std::chrono::system_clock::now();
  
  // Get Size of spacepoints
  int nSpM(0);
  int nSpB(0);
  int nSpT(0);

  for (auto sp: middleSPs) nSpM++;
  for (auto sp: bottomSPs) nSpB++;
  for (auto sp: topSPs)    nSpT++;
  if (nSpM == 0 || nSpB == 0 || nSpT == 0) return outputVec;

  CudaScalar<int> nSpM_cuda(&nSpM);
  CudaScalar<int> nSpB_cuda(&nSpB);
  CudaScalar<int> nSpT_cuda(&nSpT);

  CpuMatrix<float> spMmat_cpu(nSpM, 6); // x y z r varR varZ
  CpuMatrix<float> spBmat_cpu(nSpB, 6);
  CpuMatrix<float> spTmat_cpu(nSpT, 6);
  
  // Define Matrix and Do flattening
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > middleSPvec;
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > bottomSPvec;
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > topSPvec;  
  
  int mIdx(0);
  for (auto sp: middleSPs){
    spMmat_cpu.Set(mIdx,0,sp->x());
    spMmat_cpu.Set(mIdx,1,sp->y());
    spMmat_cpu.Set(mIdx,2,sp->z());
    spMmat_cpu.Set(mIdx,3,sp->radius());
    spMmat_cpu.Set(mIdx,4,sp->varianceR());
    spMmat_cpu.Set(mIdx,5,sp->varianceZ());
    middleSPvec.push_back(sp);
    mIdx++;
  }
  
  int bIdx(0);
  for (auto sp: bottomSPs){
    spBmat_cpu.Set(bIdx,0,sp->x());
    spBmat_cpu.Set(bIdx,1,sp->y());
    spBmat_cpu.Set(bIdx,2,sp->z());
    spBmat_cpu.Set(bIdx,3,sp->radius());
    spBmat_cpu.Set(bIdx,4,sp->varianceR());
    spBmat_cpu.Set(bIdx,5,sp->varianceZ());
    bottomSPvec.push_back(sp);
    bIdx++;
  }

  int tIdx(0);
  for (auto sp: topSPs){
    spTmat_cpu.Set(tIdx,0,sp->x());
    spTmat_cpu.Set(tIdx,1,sp->y());
    spTmat_cpu.Set(tIdx,2,sp->z());
    spTmat_cpu.Set(tIdx,3,sp->radius());
    spTmat_cpu.Set(tIdx,4,sp->varianceR());
    spTmat_cpu.Set(tIdx,5,sp->varianceZ());
    topSPvec.push_back(sp);
    tIdx++;    
  }

  CudaMatrix<float> spMmat_cuda(nSpM, 6, &spMmat_cpu); 
  CudaMatrix<float> spBmat_cuda(nSpB, 6, &spBmat_cpu);
  CudaMatrix<float> spTmat_cuda(nSpT, 6, &spTmat_cpu);
  
  /*------------------------------------
     Algorithm 1. Doublet Search (DS)
  ------------------------------------*/
  
  CudaScalar<int>  nSpMcomp_cuda(new int(0));             
  CudaScalar<int>  nSpBcompPerSpM_Max_cuda(new int(0));   
  CudaScalar<int>  nSpTcompPerSpM_Max_cuda(new int(0));   
  CudaVector<int>  nSpBcompPerSpM_cuda(nSpM); nSpBcompPerSpM_cuda.Zeros();
  CudaVector<int>  nSpTcompPerSpM_cuda(nSpM); nSpTcompPerSpM_cuda.Zeros();
  CudaVector<int>  McompIndex_cuda(nSpM);
  CudaMatrix<int>  BcompIndex_cuda(nSpB, nSpM);
  CudaMatrix<int>  TcompIndex_cuda(nSpT, nSpM);

  CudaMatrix<int>  tmpBcompIndex_cuda(nSpB, nSpM);
  CudaMatrix<int>  tmpTcompIndex_cuda(nSpT, nSpM);

  
  //dim3 DS_BlockSize;
  dim3 DS_BlockSize = m_config.maxBlockSize;
  dim3 DS_GridSize(nSpM,1,1);

  SeedfinderCudaKernels::searchDoublet( DS_GridSize, DS_BlockSize,
					nSpM_cuda.Get(), spMmat_cuda.Get(),
					nSpB_cuda.Get(), spBmat_cuda.Get(),
					nSpT_cuda.Get(), spTmat_cuda.Get(),	
					deltaRMin_cuda.Get(),
					deltaRMax_cuda.Get(), 
					cotThetaMax_cuda.Get(),
					collisionRegionMin_cuda.Get(),
					collisionRegionMax_cuda.Get(),
					nSpMcomp_cuda.Get(),
					nSpBcompPerSpM_Max_cuda.Get(),
					nSpTcompPerSpM_Max_cuda.Get(),			      
					nSpBcompPerSpM_cuda.Get(),
					nSpTcompPerSpM_cuda.Get(),
					McompIndex_cuda.Get(),
					BcompIndex_cuda.Get(),
					tmpBcompIndex_cuda.Get(),
					TcompIndex_cuda.Get(),
					tmpTcompIndex_cuda.Get()
					);  

  CpuScalar<int> nSpMcomp_cpu(&nSpMcomp_cuda);
  CpuScalar<int> nSpBcompPerSpM_Max_cpu(&nSpBcompPerSpM_Max_cuda);
  CpuScalar<int> nSpTcompPerSpM_Max_cpu(&nSpTcompPerSpM_Max_cuda);
  CpuVector<int> nSpBcompPerSpM_cpu(nSpM, &nSpBcompPerSpM_cuda);
  CpuVector<int> nSpTcompPerSpM_cpu(nSpM, &nSpTcompPerSpM_cuda);
  CpuVector<int> McompIndex_cpu(nSpM, &McompIndex_cuda);
  CpuMatrix<int> BcompIndex_cpu(nSpB, nSpM, &BcompIndex_cuda);
  CpuMatrix<int> TcompIndex_cpu(nSpT, nSpM, &TcompIndex_cuda);

  auto end_DS = std::chrono::system_clock::now();
  std::chrono::duration<double> elapse_DS = end_DS-start_DS;
  std::get<0>(t_metric) += elapse_DS.count();

  /* -----------------------------------------
     Algorithm 2. Transform coordinate
  -------------------------------------------*/
  
  auto start_TC = std::chrono::system_clock::now();  
  
  CudaMatrix<float> spMcompMat_cuda(*nSpMcomp_cpu.Get(), 6);
  
  CudaMatrix<float> spBcompMatPerSpM_cuda  (*nSpBcompPerSpM_Max_cpu.Get(), (*nSpMcomp_cpu.Get())*6);
  CudaMatrix<float> spTcompMatPerSpM_cuda  (*nSpTcompPerSpM_Max_cpu.Get(), (*nSpMcomp_cpu.Get())*6);
  CudaMatrix<float> circBcompMatPerSpM_cuda(*nSpBcompPerSpM_Max_cpu.Get(), (*nSpMcomp_cpu.Get())*6);  
  CudaMatrix<float> circTcompMatPerSpM_cuda(*nSpTcompPerSpM_Max_cpu.Get(), (*nSpMcomp_cpu.Get())*6);
  
  dim3 TC_GridSize(*nSpMcomp_cpu.Get(),1,1);
  dim3 TC_BlockSize = m_config.maxBlockSize;

  SeedfinderCudaKernels::reduceMatrix(TC_GridSize, TC_BlockSize,
				      nSpM_cuda.Get(),
				      spMmat_cuda.Get(),
				      McompIndex_cuda.Get(),
				      nSpB_cuda.Get(),
				      spBmat_cuda.Get(),
				      nSpBcompPerSpM_Max_cuda.Get(),
				      BcompIndex_cuda.Get(),
				      nSpT_cuda.Get(),
				      spTmat_cuda.Get(),
				      nSpTcompPerSpM_Max_cuda.Get(),
				      TcompIndex_cuda.Get(),
				      spMcompMat_cuda.Get(),
				      spBcompMatPerSpM_cuda.Get(),
				      circBcompMatPerSpM_cuda.Get(),
				      spTcompMatPerSpM_cuda.Get(),
				      circTcompMatPerSpM_cuda.Get()
				      );
  
  auto end_TC = std::chrono::system_clock::now();
  std::chrono::duration<double> elapse_TC = end_TC-start_TC;
  std::get<1>(t_metric) += elapse_TC.count();
  
  // ----------------------------------------------------
  //  Algorithm 3. Triplet Search (TS) & Seed filtering 
  //-----------------------------------------------------

  auto start_TS = std::chrono::system_clock::now();  
  
  // retreive middle-bottom doublet circ information
  CpuMatrix<float> circBcompMatPerSpM_cpu(*nSpBcompPerSpM_Max_cpu.Get(),
					  (*nSpMcomp_cpu.Get())*6,
					  &circBcompMatPerSpM_cuda);    

  std::vector<int> offsetVec(m_config.offsetVecSize);
  std::iota (std::begin(offsetVec), std::end(offsetVec), 0); // Fill with 0, 1, ..., 99.
  for (auto& el: offsetVec) el = el*m_config.maxBlockSize;
  CudaVector<int>   offsetVec_cuda(offsetVec.size(),&offsetVec[0]);
  
  const int         nTrplPerSpMLimit = m_config.nTrplPerSpBLimit*(*nSpBcompPerSpM_Max_cpu.Get());  
  CudaScalar<int>   nTrplPerSpMLimit_cuda(&nTrplPerSpMLimit);
  
  CudaVector<int>   nTrplPerSpM_cuda(*nSpMcomp_cpu.Get());
  nTrplPerSpM_cuda.Zeros();  
  CudaMatrix<int>   nTrplPerSpB_cuda(*nSpBcompPerSpM_Max_cpu.Get(), *nSpMcomp_cpu.Get());
  nTrplPerSpB_cuda.Zeros();  
  CudaMatrix<int>   tPassIndex_cuda(nTrplPerSpMLimit, *nSpMcomp_cpu.Get());
  CudaMatrix<int>   bPassIndex_cuda(nTrplPerSpMLimit, *nSpMcomp_cpu.Get()); 
  CudaMatrix<float> curvatures_cuda(nTrplPerSpMLimit, *nSpMcomp_cpu.Get());
  CudaMatrix<float> impactparameters_cuda(nTrplPerSpMLimit, *nSpMcomp_cpu.Get()); 
  
  CpuVector<int>    nTrplPerSpM_cpu(*nSpMcomp_cpu.Get(),true);
  nTrplPerSpM_cpu.Zeros();
  CpuMatrix<int>    nTrplPerSpB_cpu(*nSpBcompPerSpM_Max_cpu.Get(), *nSpMcomp_cpu.Get(), true);
  CpuMatrix<int>    tPassIndex_cpu(nTrplPerSpMLimit, *nSpMcomp_cpu.Get(), true);
  CpuMatrix<int>    bPassIndex_cpu(nTrplPerSpMLimit, *nSpMcomp_cpu.Get(), true); 
  CpuMatrix<float>  curvatures_cpu(nTrplPerSpMLimit, *nSpMcomp_cpu.Get(), true);
  CpuMatrix<float>  impactparameters_cpu(nTrplPerSpMLimit, *nSpMcomp_cpu.Get(),true);
  
  cudaStream_t cuStream;
  cudaStreamCreate(&cuStream);
  
  for (int i_m=0; i_m<=*nSpMcomp_cpu.Get(); i_m++){    
    
    cudaStreamSynchronize(cuStream);
    
    std::vector<std::pair<
      float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>> seedsPerSpM;

    // triplet search    
    if (i_m < *nSpMcomp_cpu.Get()){
          
      int mIndex = *McompIndex_cpu.Get(i_m);
      int nSpBcompPerSpM = *nSpBcompPerSpM_cpu.Get(mIndex);
      int nSpTcompPerSpM = *nSpTcompPerSpM_cpu.Get(mIndex);
      
      dim3 TS_GridSize(nSpBcompPerSpM,1,1);
      dim3 TS_BlockSize;
      
      int i_ts(0);    
      while ( offsetVec[i_ts] < nSpTcompPerSpM ){	
	TS_BlockSize = dim3(fmin(m_config.maxBlockSize, nSpTcompPerSpM - offsetVec[i_ts] ), 1,1); 
	SeedfinderCudaKernels::searchTriplet(TS_GridSize, TS_BlockSize,
					     offsetVec_cuda.Get(i_ts),
					     nSpMcomp_cuda.Get(),
					     spMcompMat_cuda.Get(i_m,0),
					     nSpBcompPerSpM_Max_cuda.Get(),
					     circBcompMatPerSpM_cuda.Get(0,6*i_m),
					     nSpTcompPerSpM_Max_cuda.Get(),
					     circTcompMatPerSpM_cuda.Get(offsetVec[i_ts],6*i_m),
					     // Seed finder config
					     maxScatteringAngle2_cuda.Get(),
					     sigmaScattering_cuda.Get(),
					     minHelixDiameter2_cuda.Get(),
					     pT2perRadius_cuda.Get(),
					     impactMax_cuda.Get(),
					     nTrplPerSpMLimit_cuda.Get(),
					     // output
					     nTrplPerSpM_cuda.Get(i_m),
					     nTrplPerSpB_cuda.Get(0,i_m),
					     tPassIndex_cuda.Get(0,i_m),
					     bPassIndex_cuda.Get(0,i_m),
					     curvatures_cuda.Get(0,i_m),
					     impactparameters_cuda.Get(0,i_m),
					     &cuStream
					     );
	i_ts++;
      }
      
      nTrplPerSpM_cpu.CopyD2H(nTrplPerSpM_cuda.Get(i_m),1,i_m, &cuStream);
      
      nTrplPerSpB_cpu.CopyD2H(nTrplPerSpB_cuda.Get(0,i_m),
			      *nSpBcompPerSpM_Max_cpu.Get(),
			      (*nSpBcompPerSpM_Max_cpu.Get())*i_m,
			      &cuStream);
      tPassIndex_cpu.CopyD2H(tPassIndex_cuda.Get(0,i_m),
			     nTrplPerSpMLimit,
			     nTrplPerSpMLimit*i_m,
			     &cuStream);

      bPassIndex_cpu.CopyD2H(bPassIndex_cuda.Get(0,i_m),
			     nTrplPerSpMLimit,
			     nTrplPerSpMLimit*i_m,
			     &cuStream);
      
      curvatures_cpu.CopyD2H(curvatures_cuda.Get(0,i_m),
			     nTrplPerSpMLimit,
			     nTrplPerSpMLimit*i_m,
			     &cuStream);
      
      impactparameters_cpu.CopyD2H(impactparameters_cuda.Get(0,i_m),
				   nTrplPerSpMLimit,
				   nTrplPerSpMLimit*i_m,
				   &cuStream);

    }
    
    // seed filtering
    if (i_m > 0){
      seedsPerSpM.clear();

      int curID = i_m-1;
      auto mIndex    = *McompIndex_cpu.Get(curID);
      auto bIndexVec = BcompIndex_cpu.Get(0,curID);
      auto tIndexVec = TcompIndex_cpu.Get(0,curID);
      
      // Sort Index about bottom Index
      std::map<int, std::vector< std::tuple< int, float, float > > > trplMap;

      for (int i_trpl=0; i_trpl<*nTrplPerSpM_cpu.Get(curID); i_trpl++){
	int trplBindex = *(bPassIndex_cpu.Get(i_trpl,curID));
	int trplTindex = *(tPassIndex_cpu.Get(i_trpl,curID));
	int bIndex = bIndexVec[trplBindex];
	int tIndex = tIndexVec[trplTindex];
	
	float curv     = *(curvatures_cpu.Get(i_trpl,curID));
	float impact   = *(impactparameters_cpu.Get(i_trpl,curID));

	trplMap[bIndex].push_back(std::make_tuple(tIndex, curv, impact));
      }

      std::vector<const InternalSpacePoint<external_spacepoint_t> *> tVec;
      std::vector<float> curvatures;
      std::vector<float> impactParameters;
      
      // iterate over spB
      for (auto &el: trplMap){
	tVec.clear();
	curvatures.clear();
	impactParameters.clear();

	int bIndex = el.first;
	// sort about top index
	auto triplets = el.second;
	sort(triplets.begin(), triplets.end()); 	

	for (auto it = std::make_move_iterator(triplets.begin()),
	       end = std::make_move_iterator(triplets.end()); it != end; ++it){
	  int tIndex = std::get<0>(*it);
	  tVec.push_back(std::move(topSPvec[tIndex]));	  
	  curvatures.push_back(std::move(std::get<1>(*it)));
	  impactParameters.push_back(std::move(std::get<2>(*it)));
	}

	float Zob = *(circBcompMatPerSpM_cpu.Get(el.first,(curID)*6));
	
	std::vector<std::pair<float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>> sameTrackSeeds;
	sameTrackSeeds=std::move(m_config.seedFilter->filterSeeds_2SpFixed(
								     *bottomSPvec[bIndex],
								     *middleSPvec[mIndex],
								     tVec,
								     curvatures,
								     impactParameters,Zob)); 
	seedsPerSpM.insert(seedsPerSpM.end(),
			   std::make_move_iterator(sameTrackSeeds.begin()),
			   std::make_move_iterator(sameTrackSeeds.end()));	      
	
      }
      m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, outputVec);      
    }        
  }  

  auto end_TS = std::chrono::system_clock::now();
  std::chrono::duration<double> elapse_TS = end_TS-start_TS;
  std::get<2>(t_metric) += elapse_TS.count();

  auto end_wall = std::chrono::system_clock::now();
  std::chrono::duration<double> elapse_wall = end_wall-start_wall;
  std::get<3>(t_metric) += elapse_wall.count();
  
  return outputVec;  
  }

#endif
  
}// namespace Acts

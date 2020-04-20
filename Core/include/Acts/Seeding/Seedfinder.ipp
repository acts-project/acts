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
#include <future>

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
  
  //---------------------------------
  // Algorithm 0. Matrix Flattening 
  //---------------------------------
  
  auto start_wall = std::chrono::system_clock::now();
  auto start_DS = std::chrono::system_clock::now();
  
  // Get the size of spacepoints
  int nSpM(0);
  int nSpB(0);
  int nSpT(0);

  std::for_each(middleSPs.begin(), middleSPs.end(), [&nSpM](auto sp){nSpM++;});
  std::for_each(bottomSPs.begin(), bottomSPs.end(), [&nSpB](auto sp){nSpB++;});
  std::for_each(topSPs.begin(),    topSPs.end(),    [&nSpT](auto sp){nSpT++;});

  CudaScalar<int> nSpM_cuda(&nSpM);
  CudaScalar<int> nSpB_cuda(&nSpB);
  CudaScalar<int> nSpT_cuda(&nSpT);
  
  if (nSpM == 0 || nSpB == 0 || nSpT == 0){
    return outputVec;
  }

  // Matrix flattening
  CpuMatrix<float> spMmat_cpu(nSpM, 6); // x y z r varR varZ
  CpuMatrix<float> spBmat_cpu(nSpB, 6);
  CpuMatrix<float> spTmat_cpu(nSpT, 6);
  
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > middleSPvec;
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > bottomSPvec;
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > topSPvec;  
    
  auto fillSP = [](auto& mat, auto& vec, auto& id, auto sp){
		  mat.Set(id,0,sp->x());
		  mat.Set(id,1,sp->y());
		  mat.Set(id,2,sp->z());
		  mat.Set(id,3,sp->radius());
		  mat.Set(id,4,sp->varianceR());
		  mat.Set(id,5,sp->varianceZ());
		  vec.push_back(sp);
		  id++;
		};
     
  int mIdx(0);
  for (auto sp: middleSPs){
    fillSP(spMmat_cpu, middleSPvec, mIdx, sp );
  }
  
  int bIdx(0);
  for (auto sp: bottomSPs){
    fillSP(spBmat_cpu, bottomSPvec, bIdx,sp );
  }

  int tIdx(0);
  for (auto sp: topSPs){
    fillSP(spTmat_cpu, topSPvec,   tIdx,sp );
  }

  CudaMatrix<float> spMmat_cuda(nSpM, 6, &spMmat_cpu); 
  CudaMatrix<float> spBmat_cuda(nSpB, 6, &spBmat_cpu);
  CudaMatrix<float> spTmat_cuda(nSpT, 6, &spTmat_cpu);
  
  //------------------------------------
  //  Algorithm 1. Doublet Search (DS)
  //------------------------------------
  
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

  //--------------------------------------
  //  Algorithm 2. Transform coordinate
  //--------------------------------------
  
  auto start_TC = std::chrono::system_clock::now();  
  
  CudaMatrix<float> spMcompMat_cuda(*nSpMcomp_cpu.Get(), 6);
  
  CudaMatrix<float> spBcompMatPerSpM_cuda  (*nSpBcompPerSpM_Max_cpu.Get(), (*nSpMcomp_cpu.Get())*6);
  CudaMatrix<float> spTcompMatPerSpM_cuda  (*nSpTcompPerSpM_Max_cpu.Get(), (*nSpMcomp_cpu.Get())*6);
  CudaMatrix<float> circBcompMatPerSpM_cuda(*nSpBcompPerSpM_Max_cpu.Get(), (*nSpMcomp_cpu.Get())*6);  
  CudaMatrix<float> circTcompMatPerSpM_cuda(*nSpTcompPerSpM_Max_cpu.Get(), (*nSpMcomp_cpu.Get())*6);
  
  dim3 TC_GridSize(*nSpMcomp_cpu.Get(),1,1);
  dim3 TC_BlockSize = m_config.maxBlockSize;

  SeedfinderCudaKernels::transformCoordinate(TC_GridSize, TC_BlockSize,
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
  
  //------------------------------------------------------
  //  Algorithm 3. Triplet Search (TS) & Seed filtering 
  //------------------------------------------------------

  auto start_TS = std::chrono::system_clock::now();  
  
  // retreive middle-bottom doublet circ information
  CpuMatrix<float> circBcompMatPerSpM_cpu(*nSpBcompPerSpM_Max_cpu.Get(),
					  (*nSpMcomp_cpu.Get())*6,
					  &circBcompMatPerSpM_cuda);    

  std::vector<int>  offsetVec(m_config.offsetVecSize);
  std::iota (std::begin(offsetVec), std::end(offsetVec), 0); // Fill with 0, 1, ..., 99.
  for (auto& el: offsetVec) el = el*m_config.maxBlockSize;
  CudaVector<int>   offsetVec_cuda(offsetVec.size(),&offsetVec[0]);
  
  const int         nTrplPerSpMLimit = m_config.nTrplPerSpBLimit*(*nSpBcompPerSpM_Max_cpu.Get());  
  CudaScalar<int>   nTrplPerSpMLimit_cuda(&nTrplPerSpMLimit);
  
  CudaVector<int>   nTrplPerSpM_cuda(*nSpMcomp_cpu.Get()); nTrplPerSpM_cuda.Zeros();  
  CudaMatrix<int>   TtrplIndex_cuda(nTrplPerSpMLimit, *nSpMcomp_cpu.Get());
  CudaMatrix<int>   BtrplIndex_cuda(nTrplPerSpMLimit, *nSpMcomp_cpu.Get()); 
  CudaMatrix<float> curvatures_cuda(nTrplPerSpMLimit, *nSpMcomp_cpu.Get());
  CudaMatrix<float> impactparameters_cuda(nTrplPerSpMLimit, *nSpMcomp_cpu.Get()); 
  
  CpuVector<int>    nTrplPerSpM_cpu(*nSpMcomp_cpu.Get(),true); nTrplPerSpM_cpu.Zeros();
  CpuMatrix<int>    TtrplIndex_cpu(nTrplPerSpMLimit, *nSpMcomp_cpu.Get(), true);
  CpuMatrix<int>    BtrplIndex_cpu(nTrplPerSpMLimit, *nSpMcomp_cpu.Get(), true); 
  CpuMatrix<float>  curvatures_cpu(nTrplPerSpMLimit, *nSpMcomp_cpu.Get(), true);
  CpuMatrix<float>  impactparameters_cpu(nTrplPerSpMLimit, *nSpMcomp_cpu.Get(),true);
  
  cudaStream_t cuStream;
  cudaStreamCreate(&cuStream);
  
  for (int i_m=0; i_m<=*nSpMcomp_cpu.Get(); i_m++){    
    
    cudaStreamSynchronize(cuStream);

    // Search Triplet
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
					     circTcompMatPerSpM_cuda.Get(0,6*i_m),
					     // Seed finder config
					     maxScatteringAngle2_cuda.Get(),
					     sigmaScattering_cuda.Get(),
					     minHelixDiameter2_cuda.Get(),
					     pT2perRadius_cuda.Get(),
					     impactMax_cuda.Get(),
					     nTrplPerSpMLimit_cuda.Get(),
					     // output
					     nTrplPerSpM_cuda.Get(i_m),
					     TtrplIndex_cuda.Get(0,i_m),
					     BtrplIndex_cuda.Get(0,i_m),
					     curvatures_cuda.Get(0,i_m),
					     impactparameters_cuda.Get(0,i_m),
					     &cuStream
					     );
	i_ts++;
      }
      
      nTrplPerSpM_cpu.CopyD2H(nTrplPerSpM_cuda.Get(i_m),1,i_m, &cuStream);
      
      TtrplIndex_cpu.CopyD2H(TtrplIndex_cuda.Get(0,i_m),
			     nTrplPerSpMLimit,
			     nTrplPerSpMLimit*i_m,
			     &cuStream);

      BtrplIndex_cpu.CopyD2H(BtrplIndex_cuda.Get(0,i_m),
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

    std::vector<std::pair<
      float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>> seedsPerSpM;
    
    // Seed filtering
    if (i_m > 0){

      int i_prev = i_m-1;
      auto mIndex    = *McompIndex_cpu.Get(i_prev);
      auto bIndexVec = BcompIndex_cpu.Get(0,i_prev);
      auto tIndexVec = TcompIndex_cpu.Get(0,i_prev);
      
      // Bottom index (key) is sorted in the map automatrically
      std::map<int, std::vector< std::tuple< int, float, float > > > trplMap;

      for (int i_trpl=0; i_trpl<*nTrplPerSpM_cpu.Get(i_prev); i_trpl++){
	int BtrplIndex = *(BtrplIndex_cpu.Get(i_trpl,i_prev));
	int TriplIndex = *(TtrplIndex_cpu.Get(i_trpl,i_prev));
	int bIndex = bIndexVec[BtrplIndex];
	int tIndex = tIndexVec[TriplIndex];
	
	float curv     = *(curvatures_cpu.Get(i_trpl,i_prev));
	float impact   = *(impactparameters_cpu.Get(i_trpl,i_prev));

	trplMap[bIndex].push_back(std::make_tuple(tIndex, curv, impact));
      }
      // Sort top index
      for (auto &el: trplMap){
	sort(el.second.begin(), el.second.end()); 	
      }
			            
      std::vector<const InternalSpacePoint<external_spacepoint_t> *> tVec;
      std::vector<float> curvatures;
      std::vector<float> impactParameters;
      
      for (auto &el: trplMap){
	tVec.clear();
	curvatures.clear();
	impactParameters.clear();

	int bIndex = el.first;
	auto triplets = el.second;

	for (auto it = std::make_move_iterator(triplets.begin()),
	       end = std::make_move_iterator(triplets.end()); it != end; ++it){
	  int tIndex = std::get<0>(*it);
	  tVec.push_back(std::move(topSPvec[tIndex]));	  
	  curvatures.push_back(std::move(std::get<1>(*it)));
	  impactParameters.push_back(std::move(std::get<2>(*it)));
	}

	float Zob = *(circBcompMatPerSpM_cpu.Get(bIndex,(i_prev)*6));

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

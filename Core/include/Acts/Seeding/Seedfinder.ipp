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
  }

  
  template< typename external_spacepoint_t, typename platform_t>
  template< typename T, typename sp_range_t>
  typename std::enable_if< std::is_same<T, Acts::CPU>::value, std::vector<Seed<external_spacepoint_t> > >::type
  Seedfinder<external_spacepoint_t, platform_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const {
  std::vector<Seed<external_spacepoint_t>> outputVec;

  int i_m=0;
  for (auto spM : middleSPs) {    
    i_m++;
    
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

    //std::cout << i_m << "  " << compatBottomSP.size() << "  " << compatTopSP.size() << std::endl;
    
    // contains parameters required to calculate circle with linear equation
    
    // ...for bottom-middle
    std::vector<LinCircle> linCircleBottom;
    // ...for middle-top
    std::vector<LinCircle> linCircleTop;
    
    SeedfinderCpuFunctions<external_spacepoint_t,sp_range_t>::transformCoordinates(compatBottomSP, *spM, true, linCircleBottom);
    SeedfinderCpuFunctions<external_spacepoint_t,sp_range_t>::transformCoordinates(compatTopSP, *spM, false, linCircleTop);

    auto seedsPerSpM = SeedfinderCpuFunctions<external_spacepoint_t,sp_range_t>::searchTriplet(*spM, compatBottomSP, compatTopSP, linCircleBottom, linCircleTop, m_config);
    m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, outputVec);   
  }  
  
  return outputVec;
  }

#ifdef  ACTS_HAS_CUDA
  
  // CUDA seed finding
  template< typename external_spacepoint_t, typename platform_t>
  template< typename T, typename sp_range_t>
  typename std::enable_if< std::is_same<T, Acts::CUDA>::value, std::vector<Seed<external_spacepoint_t> > >::type
  Seedfinder<external_spacepoint_t, platform_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const {
  std::vector<Seed<external_spacepoint_t>> outputVec;

  unsigned char true_cpu  = true;
  unsigned char false_cpu = false;
  CudaScalar<unsigned char> true_cuda(&true_cpu);  
  CudaScalar<unsigned char> false_cuda(&false_cpu);    
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

  // Get Size of spacepoints
  int nMiddle = 0;
  int nBottom = 0;
  int nTop    = 0;

  for (auto sp: middleSPs) nMiddle++;
  for (auto sp: bottomSPs) nBottom++;
  for (auto sp: topSPs)    nTop++;

  if (nMiddle == 0 || nBottom == 0 || nTop == 0) return outputVec;

  // Define Matrix and Do flattening
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > middleSPvec;
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > bottomSPvec;
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > topSPvec;
  
  CpuMatrix<float> spMmat_cpu(nMiddle, 6); // x y z r varR varZ
  CpuMatrix<float> spBmat_cpu(nBottom, 6);
  CpuMatrix<float> spTmat_cpu(nTop   , 6);
  
  size_t mIdx=0;
  for (auto sp: middleSPs){
    spMmat_cpu.SetEl(mIdx,0,sp->x());
    spMmat_cpu.SetEl(mIdx,1,sp->y());
    spMmat_cpu.SetEl(mIdx,2,sp->z());
    spMmat_cpu.SetEl(mIdx,3,sp->radius());
    spMmat_cpu.SetEl(mIdx,4,sp->varianceR());
    spMmat_cpu.SetEl(mIdx,5,sp->varianceZ());
    middleSPvec.push_back(sp);
    mIdx++;
  }
  CudaMatrix<float> spMmat_cuda(nMiddle, 6, &spMmat_cpu);
  CudaScalar<int>   nSpM_cuda(&nMiddle);
  
  size_t bIdx=0;
  for (auto sp: bottomSPs){
    spBmat_cpu.SetEl(bIdx,0,sp->x());
    spBmat_cpu.SetEl(bIdx,1,sp->y());
    spBmat_cpu.SetEl(bIdx,2,sp->z());
    spBmat_cpu.SetEl(bIdx,3,sp->radius());
    spBmat_cpu.SetEl(bIdx,4,sp->varianceR());
    spBmat_cpu.SetEl(bIdx,5,sp->varianceZ());
    bottomSPvec.push_back(sp);
    bIdx++;
  }

  size_t tIdx=0;
  for (auto sp: topSPs){
    spTmat_cpu.SetEl(tIdx,0,sp->x());
    spTmat_cpu.SetEl(tIdx,1,sp->y());
    spTmat_cpu.SetEl(tIdx,2,sp->z());
    spTmat_cpu.SetEl(tIdx,3,sp->radius());
    spTmat_cpu.SetEl(tIdx,4,sp->varianceR());
    spTmat_cpu.SetEl(tIdx,5,sp->varianceZ());
    topSPvec.push_back(sp);
    tIdx++;    
  }

  /*------------------------------------
     Algorithm 1. Doublet Search (DS)
  ------------------------------------*/
  
  int  BlockSize;
  dim3 DS_BlockSize;
  dim3 DS_GridSize(nMiddle,1,1);

  CudaVector<float> rM_cuda(nMiddle, spMmat_cpu.Get(0,3));
  CudaVector<float> zM_cuda(nMiddle, spMmat_cpu.Get(0,2));
  CudaVector<float> rB_cuda(nBottom, spBmat_cpu.Get(0,3));    
  CudaVector<float> zB_cuda(nBottom, spBmat_cpu.Get(0,2));
  CudaScalar<int>   nBottom_cuda(&nBottom);
  CudaVector<float> rT_cuda(nTop,    spTmat_cpu.Get(0,3));    
  CudaVector<float> zT_cuda(nTop,    spTmat_cpu.Get(0,2));  
  CudaScalar<int>   nTop_cuda(&nTop);
  
  ///// For bottom space points
  CudaMatrix<unsigned char> isCompatBottomMat_cuda(nBottom, nMiddle);

  int  offsetDS;
  offsetDS=0;
  while(offsetDS<nBottom){
    DS_BlockSize = dim3(fmin(MAX_BLOCK_SIZE, nBottom-offsetDS), 1,1);
    SeedfinderCudaKernels::searchDoublet( DS_GridSize, DS_BlockSize,
					  true_cuda.Get(),
					  rM_cuda.Get(), zM_cuda.Get(),
					  nBottom_cuda.Get(), rB_cuda.Get(offsetDS), zB_cuda.Get(offsetDS), 
					  deltaRMin_cuda.Get(), deltaRMax_cuda.Get(), 
					  cotThetaMax_cuda.Get(),
					  collisionRegionMin_cuda.Get(),collisionRegionMax_cuda.Get(),
					  isCompatBottomMat_cuda.Get(offsetDS,0));
    offsetDS+=DS_BlockSize.x;
  }
  CpuMatrix<unsigned char>  isCompatBottomMat_cpu(nBottom, nMiddle, &isCompatBottomMat_cuda);
  
  ///// For top space points
  CudaMatrix<unsigned char> isCompatTopMat_cuda(nTop, nMiddle);
  offsetDS=0;
  while(offsetDS<nTop){
    DS_BlockSize = dim3(fmin(MAX_BLOCK_SIZE, nTop-offsetDS), 1,1);    
    SeedfinderCudaKernels::searchDoublet( DS_GridSize, DS_BlockSize,
					  false_cuda.Get(),
					  rM_cuda.Get(), zM_cuda.Get(),
					  nTop_cuda.Get(), rT_cuda.Get(offsetDS), zT_cuda.Get(offsetDS), 
					  deltaRMin_cuda.Get(), deltaRMax_cuda.Get(), 
					  cotThetaMax_cuda.Get(),
					  collisionRegionMin_cuda.Get(),collisionRegionMax_cuda.Get(),
					  isCompatTopMat_cuda.Get(offsetDS,0));
    offsetDS+=DS_BlockSize.x;
  }
  CpuMatrix<unsigned char>  isCompatTopMat_cpu(nTop, nMiddle, &isCompatTopMat_cuda);

  /* -----------------------------------------
     Algorithm 2. Transform Coordinates (TC)
  -------------------------------------------*/
  
  std::vector< std::tuple<int, std::vector< int >, std::vector< int > > > mCompIndex;

  int nBcompMax = 0;
  int nTcompMax = 0;

  for (int i_m=0; i_m<nMiddle; i_m++){
    std::vector< int > bIndex;
    for (int i=0; i<nBottom; i++){
      if (*isCompatBottomMat_cpu.Get(i,i_m)) bIndex.push_back(i);
    }
    if (bIndex.empty()) continue;    
    std::vector< int > tIndex;
    for (int i=0; i<nTop; i++){
      if (*isCompatTopMat_cpu.Get(i,i_m)) tIndex.push_back(i);
    }
    if (tIndex.empty()) continue;
    
    auto tup = std::make_tuple(i_m, bIndex, tIndex);
    mCompIndex.push_back(tup);
    nBcompMax = fmax(bIndex.size(), nBcompMax);
    nTcompMax = fmax(tIndex.size(), nTcompMax);
  }
  
  CpuMatrix<float>  spMcompMat_cpu(mCompIndex.size(), 6);
  CudaMatrix<float> spMcompMat_cuda(mCompIndex.size(),6);
  CudaScalar<int>   nBcompMax_cuda(&nBcompMax);
  CudaScalar<int>   nTcompMax_cuda(&nTcompMax);

  CpuMatrix<float>  spBcompMat_cpu   (nBcompMax, mCompIndex.size()*6);
  CudaMatrix<float> spBcompMat_cuda  (nBcompMax, mCompIndex.size()*6);
  CudaMatrix<float> circBcompMat_cuda(nBcompMax, mCompIndex.size()*6);  
  CpuMatrix<float>  spTcompMat_cpu   (nTcompMax, mCompIndex.size()*6);  
  CudaMatrix<float> spTcompMat_cuda  (nTcompMax, mCompIndex.size()*6);   
  CudaMatrix<float> circTcompMat_cuda(nTcompMax, mCompIndex.size()*6);
  
  for (int i_m=0; i_m<spMcompMat_cpu.GetNRows(); i_m++){
    auto mIndex = std::get<0>(mCompIndex[i_m]);
    auto bIndex = std::get<1>(mCompIndex[i_m]);
    auto tIndex = std::get<2>(mCompIndex[i_m]);

    //std::cout << mIndex << "  " << bIndex.size() << "  " << tIndex.size() << std::endl;
    
    for (int i=0; i<6; i++){
      spMcompMat_cpu.SetEl(i_m,i,std::move(*spMmat_cpu.Get(mIndex,i)));
    }      
    //spMcompMat_cpu.SetEl(i_m,0,*spMmat_cpu.Get(mIndex,0));
    //spMcompMat_cpu.SetEl(i_m,1,*spMmat_cpu.Get(mIndex,1));
    //spMcompMat_cpu.SetEl(i_m,2,*spMmat_cpu.Get(mIndex,2));
    //spMcompMat_cpu.SetEl(i_m,3,*spMmat_cpu.Get(mIndex,3));
    //spMcompMat_cpu.SetEl(i_m,4,*spMmat_cpu.Get(mIndex,4));
    //spMcompMat_cpu.SetEl(i_m,5,*spMmat_cpu.Get(mIndex,5));
    
    for (int i_b=0; i_b<bIndex.size(); i_b++){
      for (int i=0; i<6; i++){
	spBcompMat_cpu.SetEl(i_b, i_m*6+i, std::move(*spBmat_cpu.Get(bIndex[i_b],i)));
      }            
      //spBcompMat_cpu.SetEl(i_b, i_m*6+0, *spBmat_cpu.Get(bIndex[i_b],0));
      //spBcompMat_cpu.SetEl(i_b, i_m*6+1, *spBmat_cpu.Get(bIndex[i_b],1));
      //spBcompMat_cpu.SetEl(i_b, i_m*6+2, *spBmat_cpu.Get(bIndex[i_b],2));	    
      //spBcompMat_cpu.SetEl(i_b, i_m*6+3, *spBmat_cpu.Get(bIndex[i_b],3));
      //spBcompMat_cpu.SetEl(i_b, i_m*6+4, *spBmat_cpu.Get(bIndex[i_b],4));
      //spBcompMat_cpu.SetEl(i_b, i_m*6+5, *spBmat_cpu.Get(bIndex[i_b],5));
    }

    for (int i_t=0; i_t<tIndex.size(); i_t++){
      for (int i=0; i<6; i++){
	spTcompMat_cpu.SetEl(i_t, i_m*6+i, std::move(*spTmat_cpu.Get(tIndex[i_t],i)));
      }                  
      //spTcompMat_cpu.SetEl(i_t, i_m*6+0, *spTmat_cpu.Get(tIndex[i_t],0));
      //spTcompMat_cpu.SetEl(i_t, i_m*6+1, *spTmat_cpu.Get(tIndex[i_t],1));
      //spTcompMat_cpu.SetEl(i_t, i_m*6+2, *spTmat_cpu.Get(tIndex[i_t],2));	    
      //spTcompMat_cpu.SetEl(i_t, i_m*6+3, *spTmat_cpu.Get(tIndex[i_t],3));
      //spTcompMat_cpu.SetEl(i_t, i_m*6+4, *spTmat_cpu.Get(tIndex[i_t],4));
      //spTcompMat_cpu.SetEl(i_t, i_m*6+5, *spTmat_cpu.Get(tIndex[i_t],5));
    }    
  }

  spMcompMat_cuda.CopyH2D(spMcompMat_cpu.Get(), spMcompMat_cpu.GetSize(), 0);
  spBcompMat_cuda.CopyH2D(spBcompMat_cpu.Get(), spBcompMat_cpu.GetSize(), 0);
  spTcompMat_cuda.CopyH2D(spTcompMat_cpu.Get(), spTcompMat_cpu.GetSize(), 0);

  dim3 TC_GridSize(spMcompMat_cpu.GetNRows(),1,1);
  dim3 TC_BlockSize;
  
  // For bottom-middle
  int offsetTC;
  offsetTC=0;
  while(offsetTC<nBcompMax){
    TC_BlockSize = dim3(fmin(MAX_BLOCK_SIZE, nBcompMax-offsetTC),1,1);    
    SeedfinderCudaKernels::transformCoordinates(TC_GridSize, TC_BlockSize,
						true_cuda.Get(),
						spMcompMat_cuda.Get(0,0),
						nBcompMax_cuda.Get(),
						spBcompMat_cuda.Get(offsetTC,0),
						circBcompMat_cuda.Get(offsetTC,0));    
    offsetTC+=TC_BlockSize.x;
  }
  
  // For middle-top 
  offsetTC=0;
  while(offsetTC<nTcompMax){
    TC_BlockSize = dim3(fmin(MAX_BLOCK_SIZE, nTcompMax-offsetTC),1,1);
    SeedfinderCudaKernels::transformCoordinates(TC_GridSize, TC_BlockSize,
						false_cuda.Get(),
						spMcompMat_cuda.Get(0,0),
						nTcompMax_cuda.Get(),
						spTcompMat_cuda.Get(offsetTC,0),
						circTcompMat_cuda.Get(offsetTC,0));    
    offsetTC+=TC_BlockSize.x;
  }
  
  int nMcomp = mCompIndex.size();  
  CudaScalar<int>  nMcomp_cuda(&nMcomp);
  
  // retreive middle-bottom doublet circ information
  CpuMatrix<float> circBcompMat_cpu(nBcompMax, mCompIndex.size()*6);
  circBcompMat_cpu.CopyD2H(circBcompMat_cuda.Get(0,0), circBcompMat_cuda.GetSize(), 0);  

  std::vector<int> offsetVec(m_config.offsetVecSize);
  std::iota (std::begin(offsetVec), std::end(offsetVec), 0); // Fill with 0, 1, ..., 99.
  for (auto& el: offsetVec) el = el*MAX_BLOCK_SIZE;
  CudaVector<int> offsetVec_cuda(offsetVec.size(),&offsetVec[0]);
    
  const int         nTopPassLimit = m_config.nTopPassLimit;
  CudaScalar<int>   nTopPassLimit_cuda(&nTopPassLimit);
  
  CudaMatrix<int>   nTopPass_cuda(nBcompMax, nMcomp);
  CudaMatrix<int>   tPassIndex_cuda(nTopPassLimit, nBcompMax*nMcomp); 
  CudaMatrix<float> curvatures_cuda(nTopPassLimit, nBcompMax*nMcomp);       
  CudaMatrix<float> impactparameters_cuda(nTopPassLimit, nBcompMax*nMcomp);

  CpuMatrix<int>    nTopPass_cpu(nBcompMax, nMcomp);
  CpuMatrix<int>    tPassIndex_cpu(nTopPassLimit, nBcompMax*nMcomp); 
  CpuMatrix<float>  curvatures_cpu(nTopPassLimit, nBcompMax*nMcomp);
  CpuMatrix<float>  impactparameters_cpu(nTopPassLimit, nBcompMax*nMcomp);

  cudaStream_t cuStream;
  cudaStreamCreate(&cuStream);
   
  for (int i_m=0; i_m<=mCompIndex.size(); i_m++){
    
    cudaStreamSynchronize(cuStream);

    // -----------------------------------
    //  Algorithm 3. Triplet Search (TS)
    //------------------------------------
    
    std::vector<std::pair<
      float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>> seedsPerSpM;
    
    std::vector<const InternalSpacePoint<external_spacepoint_t> *> tVec;
    std::vector<float> curvatures;
    std::vector<float> impactParameters;
    
    if (i_m < mCompIndex.size()){
          
      // SEED Filtering is done ASYNCHRONOUSLY against Triplet search
      // Need to call it again after last iteration
            
      // For triplets collected at previous iteration
      
      auto mIndex = std::get<0>(mCompIndex[i_m]);
      auto bIndex = std::get<1>(mCompIndex[i_m]);
      auto tIndex = std::get<2>(mCompIndex[i_m]);
      
      dim3 TS_GridSize(bIndex.size(),1,1);
      dim3 TS_BlockSize;
      
      int i_ts = 0;    
      while ( offsetVec[i_ts] < tIndex.size() ){
	
	TS_BlockSize = dim3(fmin(MAX_BLOCK_SIZE, tIndex.size()-offsetVec[i_ts] ),
			    1,1);
      	
	SeedfinderCudaKernels::searchTriplet(TS_GridSize, TS_BlockSize,
					     offsetVec_cuda.Get(i_ts),
					     nMcomp_cuda.Get(),
					     spMcompMat_cuda.Get(i_m,0),
					     nBcompMax_cuda.Get(),
					     spBcompMat_cuda.Get(0,6*i_m),
					     nTcompMax_cuda.Get(),
					     spTcompMat_cuda.Get(offsetVec[i_ts],6*i_m),
					     circBcompMat_cuda.Get(0,6*i_m),
					     circTcompMat_cuda.Get(offsetVec[i_ts],6*i_m),
					     // Seed finder config
					     maxScatteringAngle2_cuda.Get(),
					     sigmaScattering_cuda.Get(),
					     minHelixDiameter2_cuda.Get(),
					     pT2perRadius_cuda.Get(),
					     impactMax_cuda.Get(),
					     nTopPassLimit_cuda.Get(),
					     // output
					     nTopPass_cuda.Get(0,i_m),
					     tPassIndex_cuda.Get(0,nBcompMax*i_m),
					     curvatures_cuda.Get(0,nBcompMax*i_m),
					     impactparameters_cuda.Get(0,nBcompMax*i_m),
					     &cuStream
					     );
	i_ts++;
      }
      
      nTopPass_cpu.CopyD2H(nTopPass_cuda.Get(0,i_m),
			   nBcompMax,
			   nBcompMax*i_m,
			   &cuStream);
      tPassIndex_cpu.CopyD2H(tPassIndex_cuda.Get(0,nBcompMax*i_m),
			     nTopPassLimit*nBcompMax,
			     nTopPassLimit*nBcompMax*i_m,
			     &cuStream);
      curvatures_cpu.CopyD2H(curvatures_cuda.Get(0,nBcompMax*i_m),
			     nTopPassLimit*nBcompMax,
			     nTopPassLimit*nBcompMax*i_m,
			     &cuStream); 
      impactparameters_cpu.CopyD2H(impactparameters_cuda.Get(0,nBcompMax*i_m),
				   nTopPassLimit*nBcompMax,
				   nTopPassLimit*nBcompMax*i_m,
				   &cuStream);
    }

    // --------------------------------
    //  Algorithm 4. Seed Filter (SF)
    // --------------------------------
    
    if (i_m > 0){
      seedsPerSpM.clear();
      auto middleIdx     = std::get<0>(mCompIndex[i_m-1]);
      auto compBottomIdx = std::get<1>(mCompIndex[i_m-1]);
      auto compTopIdx    = std::get<2>(mCompIndex[i_m-1]);
      
      for (int i_b=0; i_b<compBottomIdx.size(); i_b++){
	int nTpass = *(nTopPass_cpu.Get(i_b,i_m-1));

	if (nTpass==0) continue;	
	
	tVec.clear();
	curvatures.clear();
	impactParameters.clear();      
	float Zob = *(circBcompMat_cpu.Get(i_b,(i_m-1)*6));

	std::vector< std::tuple< int, int, int > > indexVec;
	for(int i_t=0; i_t<nTpass; i_t++){
	  int g_tIndex = compTopIdx[*tPassIndex_cpu.Get(i_t,i_b+(i_m-1)*nBcompMax)];
	  indexVec.push_back(std::make_tuple(g_tIndex,i_t,i_b));
	}
	sort(indexVec.begin(), indexVec.end()); 
	
	for(auto el: indexVec){
	  auto g_tIndex = std::get<0>(el);
	  auto tId      = std::get<1>(el);
	  auto bId      = std::get<2>(el);
	  
	  tVec.push_back(topSPvec[g_tIndex]);
	  curvatures.push_back(*curvatures_cpu.Get(tId,bId+(i_m-1)*nBcompMax));
	  impactParameters.push_back(*impactparameters_cpu.Get(tId,bId+(i_m-1)*nBcompMax));
	}
	
	std::vector<std::pair<float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>> sameTrackSeeds;
	sameTrackSeeds = std::move(m_config.seedFilter->filterSeeds_2SpFixed(*bottomSPvec[compBottomIdx[i_b]],
									     *middleSPvec[middleIdx],
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
  
  return outputVec;  
  }

#endif
  
}// namespace Acts

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
    spMmat_cpu.Set(mIdx,0,sp->x());
    spMmat_cpu.Set(mIdx,1,sp->y());
    spMmat_cpu.Set(mIdx,2,sp->z());
    spMmat_cpu.Set(mIdx,3,sp->radius());
    spMmat_cpu.Set(mIdx,4,sp->varianceR());
    spMmat_cpu.Set(mIdx,5,sp->varianceZ());
    middleSPvec.push_back(sp);
    mIdx++;
  }
  
  size_t bIdx=0;
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

  size_t tIdx=0;
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

  CudaMatrix<float> spMmat_cuda(nMiddle, 6, &spMmat_cpu); 
  CudaMatrix<float> spBmat_cuda(nBottom, 6, &spBmat_cpu);
  CudaMatrix<float> spTmat_cuda(nTop   , 6, &spTmat_cpu);
  CudaScalar<int>   nSpM_cuda(&nMiddle);
  CudaScalar<int>   nSpB_cuda(&nBottom);
  CudaScalar<int>   nSpT_cuda(&nTop);
  
  /*------------------------------------
     Algorithm 1. Doublet Search (DS)
  ------------------------------------*/
  
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

  CudaVector<int>   nSpBcomp_cuda(nMiddle); nSpBcomp_cuda.Zeros();
  CudaVector<int>   nSpTcomp_cuda(nMiddle); nSpTcomp_cuda.Zeros();
  
  // For bottom space points
  CudaMatrix<bool> isCompatBottomMat_cuda(nBottom, nMiddle);

  int offsetDS;
  offsetDS=0;
  while(offsetDS<nBottom){
    DS_BlockSize = dim3(fmin(m_config.maxBlockSize, nBottom-offsetDS), 1,1);
    SeedfinderCudaKernels::searchDoublet( DS_GridSize, DS_BlockSize,
					  true_cuda.Get(),
					  rM_cuda.Get(), zM_cuda.Get(),
					  nBottom_cuda.Get(),
					  rB_cuda.Get(offsetDS),
					  zB_cuda.Get(offsetDS), 
					  deltaRMin_cuda.Get(), deltaRMax_cuda.Get(), 
					  cotThetaMax_cuda.Get(),
					  collisionRegionMin_cuda.Get(),collisionRegionMax_cuda.Get(),
					  isCompatBottomMat_cuda.Get(offsetDS,0),
					  nSpBcomp_cuda.Get());
    offsetDS+=DS_BlockSize.x;
  }
  
  // For top space points
  CudaMatrix<bool> isCompatTopMat_cuda(nTop, nMiddle);
  offsetDS=0;
  while(offsetDS<nTop){
    DS_BlockSize = dim3(fmin(m_config.maxBlockSize, nTop-offsetDS), 1,1);    
    SeedfinderCudaKernels::searchDoublet( DS_GridSize, DS_BlockSize,
					  false_cuda.Get(),
					  rM_cuda.Get(), zM_cuda.Get(),
					  nTop_cuda.Get(),
					  rT_cuda.Get(offsetDS),
					  zT_cuda.Get(offsetDS), 
					  deltaRMin_cuda.Get(), deltaRMax_cuda.Get(), 
					  cotThetaMax_cuda.Get(),
					  collisionRegionMin_cuda.Get(),collisionRegionMax_cuda.Get(),
					  isCompatTopMat_cuda.Get(offsetDS,0),
					  nSpTcomp_cuda.Get());
    offsetDS+=DS_BlockSize.x;
  }

  CpuMatrix<bool>  isCompatBottomMat_cpu(nBottom, nMiddle, &isCompatBottomMat_cuda);
  CpuMatrix<bool>  isCompatTopMat_cpu(nTop, nMiddle, &isCompatTopMat_cuda);
  auto nSpBcomp_cpu = nSpBcomp_cuda.GetHost();
  auto nSpTcomp_cpu = nSpTcomp_cuda.GetHost();
  
  auto end_DS = std::chrono::system_clock::now();
  std::chrono::duration<double> elapse_DS = end_DS-start_DS;
  std::get<0>(t_metric) += elapse_DS.count();
 
  // Count number of compatible hits
  int nMcomp = 0;
  int nBcompMax = 0;
  int nTcompMax = 0;

  //int nSpMcomp = 0;
  //int nSpBcompPerSpM_Max = 0;
  //int nSpTcompPerSpM_Max = 0;
  
  for (int i_m=0; i_m<nMiddle; i_m++){
    if ( !nSpBcomp_cpu[i_m] || !nSpTcomp_cpu[i_m] ) continue;
    nBcompMax = fmax(nSpBcomp_cpu[i_m], nBcompMax);
    nTcompMax = fmax(nSpTcomp_cpu[i_m], nTcompMax);    
    nMcomp++;
  }

  CudaScalar<int>   nMcomp_cuda(&nMcomp);
  CudaScalar<int>   nBcompMax_cuda(&nBcompMax);
  CudaScalar<int>   nTcompMax_cuda(&nTcompMax);

  // Get the index of compatible hits
  auto mIndex_cpu = new int[nMcomp];  
  CpuMatrix< int > bIndex_cpu(nBcompMax, nMcomp); bIndex_cpu.Zeros();
  CpuMatrix< int > tIndex_cpu(nTcompMax, nMcomp); tIndex_cpu.Zeros();
    
  int mcnt=0;
  for (int i_m=0; i_m<nMiddle; i_m++){
    if ( !nSpBcomp_cpu[i_m] || !nSpTcomp_cpu[i_m] ) continue;

    auto BottomCol = isCompatBottomMat_cpu.Get(0,i_m);
    int bcnt = 0;
    for (int i=0; i<nBottom; i++){
      if (BottomCol[i]){
	bIndex_cpu.Set(bcnt,mcnt,i);
	bcnt++;
      }
    }
    
    auto TopCol = isCompatTopMat_cpu.Get(0,i_m);    
    int tcnt = 0;
    for (int i=0; i<nTop; i++){
      if (TopCol[i]){
	tIndex_cpu.Set(tcnt,mcnt,i);
	tcnt++;
      }
    }

    mIndex_cpu[mcnt] = i_m;    
    mcnt++;    
  }
  
  CudaVector< int > mIndex_cuda(nMcomp,    mIndex_cpu);
  CudaMatrix< int > bIndex_cuda(nBcompMax, nMcomp, &bIndex_cpu);
  CudaMatrix< int > tIndex_cuda(nTcompMax, nMcomp, &tIndex_cpu);

  /* -----------------------------------------
     Algorithm 1.B Matrix reduction
  -------------------------------------------*/
  
  CudaMatrix<float> spBcompMat_cuda(nBcompMax, nMcomp*6);
  CudaMatrix<float> spTcompMat_cuda(nTcompMax, nMcomp*6);

  dim3 RM_GridSize(nMcomp,1,1);
  dim3 RM_BlockSize; 

  // For bottom hits
  int offsetRM;
  offsetRM=0;
  while(offsetRM<nBcompMax){
    RM_BlockSize = dim3(fmin(m_config.maxBlockSize, nBcompMax-offsetRM),1,1);
    SeedfinderCudaKernels::reduceMatrix(RM_GridSize, RM_BlockSize,
					nSpB_cuda.Get(),
					spBmat_cuda.Get(0,0),
					nBcompMax_cuda.Get(),
					bIndex_cuda.Get(offsetRM,0),
					spBcompMat_cuda.Get(offsetRM,0)
					);
    offsetRM+=RM_BlockSize.x;    
  }

  // For top hits
  offsetRM=0;
  while(offsetRM<nTcompMax){
    RM_BlockSize = dim3(fmin(m_config.maxBlockSize, nTcompMax-offsetRM),1,1);
    SeedfinderCudaKernels::reduceMatrix(RM_GridSize, RM_BlockSize,
					nSpT_cuda.Get(),
					spTmat_cuda.Get(0,0),
					nTcompMax_cuda.Get(),
					tIndex_cuda.Get(offsetRM,0),
					spTcompMat_cuda.Get(offsetRM,0)
					);
    offsetRM+=RM_BlockSize.x;    
  }

  // For middle hits
  CpuMatrix<float>  spMcompMat_cpu (nMcomp, 6);
  CudaMatrix<float> spMcompMat_cuda(nMcomp, 6);
  for (int i_m=0; i_m<nMcomp; i_m++){
    int mIndex = mIndex_cpu[i_m];    
    spMcompMat_cpu.Set(i_m,0,*spMmat_cpu.Get(mIndex,0));
    spMcompMat_cpu.Set(i_m,1,*spMmat_cpu.Get(mIndex,1));
    spMcompMat_cpu.Set(i_m,2,*spMmat_cpu.Get(mIndex,2));
    spMcompMat_cpu.Set(i_m,3,*spMmat_cpu.Get(mIndex,3));
    spMcompMat_cpu.Set(i_m,4,*spMmat_cpu.Get(mIndex,4));
    spMcompMat_cpu.Set(i_m,5,*spMmat_cpu.Get(mIndex,5));    
  }
  spMcompMat_cuda.CopyH2D(spMcompMat_cpu.Get(), spMcompMat_cpu.GetSize(), 0);

  /* -----------------------------------------
     Algorithm 2. Transform Coordinates (TC)
  -------------------------------------------*/

  auto start_TC = std::chrono::system_clock::now();  
  
  CudaMatrix<float> circBcompMat_cuda(nBcompMax, nMcomp*6);  
  CudaMatrix<float> circTcompMat_cuda(nTcompMax, nMcomp*6);
    
  dim3 TC_GridSize(nMcomp,1,1);
  dim3 TC_BlockSize;
  
  // For bottom-middle doublets
  int offsetTC;
  offsetTC=0;
  while(offsetTC<nBcompMax){
    TC_BlockSize = dim3(fmin(m_config.maxBlockSize, nBcompMax-offsetTC),1,1);    
    SeedfinderCudaKernels::transformCoordinates(TC_GridSize,
						TC_BlockSize,
						true_cuda.Get(),
						spMcompMat_cuda.Get(0,0),
						nBcompMax_cuda.Get(),
						spBcompMat_cuda.Get(offsetTC,0),
						circBcompMat_cuda.Get(offsetTC,0));    
    offsetTC+=TC_BlockSize.x;
  }
  
  // For middle-top doublets
  offsetTC=0;
  while(offsetTC<nTcompMax){
    TC_BlockSize = dim3(fmin(m_config.maxBlockSize, nTcompMax-offsetTC),1,1);
    SeedfinderCudaKernels::transformCoordinates(TC_GridSize,
						TC_BlockSize,
						false_cuda.Get(),
						spMcompMat_cuda.Get(0,0),
						nTcompMax_cuda.Get(),
						spTcompMat_cuda.Get(offsetTC,0),
						circTcompMat_cuda.Get(offsetTC,0));    
    offsetTC+=TC_BlockSize.x;
  }
    
  // retreive middle-bottom doublet circ information
  CpuMatrix<float> circBcompMat_cpu(nBcompMax, nMcomp*6, &circBcompMat_cuda);
  
  auto end_TC = std::chrono::system_clock::now();
  std::chrono::duration<double> elapse_TC = end_TC-start_TC;
  std::get<1>(t_metric) += elapse_TC.count();

  auto start_TS = std::chrono::system_clock::now();  
  
  std::vector<int> offsetVec(m_config.offsetVecSize);
  std::iota (std::begin(offsetVec), std::end(offsetVec), 0); // Fill with 0, 1, ..., 99.
  for (auto& el: offsetVec) el = el*m_config.maxBlockSize;
  CudaVector<int> offsetVec_cuda(offsetVec.size(),&offsetVec[0]);
  
  const int         nTopPassLimit = m_config.nTopPassLimit;
  CudaScalar<int>   nTopPassLimit_cuda(&nTopPassLimit);

  CudaMatrix<int>   nTopPass_cuda(nBcompMax, nMcomp); nTopPass_cuda.Zeros();
  CudaMatrix<int>   tPassIndex_cuda(nTopPassLimit, nBcompMax*nMcomp); 
  CudaMatrix<float> curvatures_cuda(nTopPassLimit, nBcompMax*nMcomp);       
  CudaMatrix<float> impactparameters_cuda(nTopPassLimit, nBcompMax*nMcomp);

  CpuMatrix<int>    nTopPass_cpu(nBcompMax, nMcomp, true);
  CpuMatrix<int>    tPassIndex_cpu(nTopPassLimit, nBcompMax*nMcomp, true);
  CpuMatrix<float>  curvatures_cpu(nTopPassLimit, nBcompMax*nMcomp, true);
  CpuMatrix<float>  impactparameters_cpu(nTopPassLimit, nBcompMax*nMcomp, true);
  
  cudaStream_t cuStream;
  cudaStreamCreate(&cuStream);
   
  for (int i_m=0; i_m<=nMcomp; i_m++){    
    cudaStreamSynchronize(cuStream);
    // -----------------------------------
    //  Algorithm 3. Triplet Search (TS)
    //------------------------------------
    
    std::vector<std::pair<
      float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>> seedsPerSpM;
    
    std::vector<const InternalSpacePoint<external_spacepoint_t> *> tVec;
    std::vector<float> curvatures;
    std::vector<float> impactParameters;
    
    if (i_m < nMcomp){
          
      // For triplets collected at the previous iteration      
      int mIndex = mIndex_cpu[i_m];      
      dim3 TS_GridSize(nSpBcomp_cpu[mIndex],1,1);
      dim3 TS_BlockSize;
      
      int i_ts = 0;    
      while ( offsetVec[i_ts] < nSpTcomp_cpu[mIndex] ){	
	TS_BlockSize = dim3(fmin(m_config.maxBlockSize, nSpTcomp_cpu[mIndex]-offsetVec[i_ts] ), 1,1); 
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
      int curID = i_m-1;
      seedsPerSpM.clear();

      int mIndex = mIndex_cpu[curID]; 
      auto bIndexVec = bIndex_cpu.Get(0,curID);
      auto tIndexVec = tIndex_cpu.Get(0,curID);
      
      for (int i_b=0; i_b<nSpBcomp_cpu[mIndex]; i_b++){
	int nTpass = *(nTopPass_cpu.Get(i_b,curID));
	if (nTpass==0) continue;	
	
	tVec.clear();
	curvatures.clear();
	impactParameters.clear();      
	float Zob = *(circBcompMat_cpu.Get(i_b,(curID)*6));

	std::vector< std::tuple< int, int, int > > indexVec;
	for(int i_t=0; i_t<nTpass; i_t++){
	  int g_tIndex = tIndexVec[*tPassIndex_cpu.Get(i_t,i_b+(curID)*nBcompMax)];
	  indexVec.push_back(std::make_tuple(g_tIndex,i_t,i_b));
	}
	sort(indexVec.begin(), indexVec.end()); 
	
	for(auto el: indexVec){
	  auto g_tIndex = std::get<0>(el);
	  auto tId      = std::get<1>(el);
	  auto bId      = std::get<2>(el);
	  
	  tVec.push_back(topSPvec[g_tIndex]);
	  curvatures.push_back(*curvatures_cpu.Get(tId,bId+(curID)*nBcompMax));
	  impactParameters.push_back(*impactparameters_cpu.Get(tId,bId+(curID)*nBcompMax));
	}
	
	std::vector<std::pair<float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>> sameTrackSeeds;
	sameTrackSeeds = std::move(m_config.seedFilter->filterSeeds_2SpFixed(*bottomSPvec[bIndexVec[i_b]],
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

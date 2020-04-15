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
  
  // Define Matrix and Do flattening
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > middleSPvec;
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > bottomSPvec;
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > topSPvec;
  
  CpuMatrix<float> spMmat_cpu(nSpM, 6); // x y z r varR varZ
  CpuMatrix<float> spBmat_cpu(nSpB, 6);
  CpuMatrix<float> spTmat_cpu(nSpT, 6);
  
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
  
  dim3 DS_BlockSize;
  dim3 DS_GridSize(nSpM,1,1);

  CudaVector<int>   nSpBcompPerSpM_cuda(nSpM); nSpBcompPerSpM_cuda.Zeros();
  CudaVector<int>   nSpTcompPerSpM_cuda(nSpM); nSpTcompPerSpM_cuda.Zeros();
  
  // For bottom space points
  CudaMatrix<bool> isCompatBottomMat_cuda(nSpB, nSpM);

  int offsetDS;
  offsetDS=0;
  while(offsetDS<nSpB){
    DS_BlockSize = dim3(fmin(m_config.maxBlockSize, nSpB-offsetDS), 1,1);
    SeedfinderCudaKernels::searchDoublet( DS_GridSize, DS_BlockSize,
					  true_cuda.Get(),
					  nSpM_cuda.Get(), spMmat_cuda.Get(),
					  nSpB_cuda.Get(), spBmat_cuda.Get(offsetDS,0),
					  deltaRMin_cuda.Get(), deltaRMax_cuda.Get(), 
					  cotThetaMax_cuda.Get(),
					  collisionRegionMin_cuda.Get(),collisionRegionMax_cuda.Get(),
					  isCompatBottomMat_cuda.Get(offsetDS,0),
					  nSpBcompPerSpM_cuda.Get());
    offsetDS+=DS_BlockSize.x;
  }
  
  // For top space points
  CudaMatrix<bool> isCompatTopMat_cuda(nSpT, nSpM);
  offsetDS=0;
  while(offsetDS<nSpT){
    DS_BlockSize = dim3(fmin(m_config.maxBlockSize, nSpT-offsetDS), 1,1);    
    SeedfinderCudaKernels::searchDoublet( DS_GridSize, DS_BlockSize,
					  false_cuda.Get(),
					  nSpM_cuda.Get(), spMmat_cuda.Get(),
					  nSpT_cuda.Get(), spTmat_cuda.Get(offsetDS,0),
					  deltaRMin_cuda.Get(), deltaRMax_cuda.Get(), 
					  cotThetaMax_cuda.Get(),
					  collisionRegionMin_cuda.Get(),collisionRegionMax_cuda.Get(),
					  isCompatTopMat_cuda.Get(offsetDS,0),
					  nSpTcompPerSpM_cuda.Get());
    offsetDS+=DS_BlockSize.x;
  }

  CpuMatrix<bool>  isCompatBottomMat_cpu(nSpB, nSpM, &isCompatBottomMat_cuda);
  CpuMatrix<bool>  isCompatTopMat_cpu(nSpT, nSpM, &isCompatTopMat_cuda);
  auto nSpBcompPerSpM_cpu = nSpBcompPerSpM_cuda.GetHost();
  auto nSpTcompPerSpM_cpu = nSpTcompPerSpM_cuda.GetHost();
  
  auto end_DS = std::chrono::system_clock::now();
  std::chrono::duration<double> elapse_DS = end_DS-start_DS;
  std::get<0>(t_metric) += elapse_DS.count();
 
  // Count number of compatible hits
  int nSpMcomp(0);
  int nSpBcompPerSpM_Max(0);
  int nSpTcompPerSpM_Max(0);
  
  for (int i_m=0; i_m<nSpM; i_m++){
    if ( !nSpBcompPerSpM_cpu[i_m] || !nSpTcompPerSpM_cpu[i_m] ) continue;

    nSpBcompPerSpM_Max = fmax(nSpBcompPerSpM_cpu[i_m], nSpBcompPerSpM_Max);
    nSpTcompPerSpM_Max = fmax(nSpTcompPerSpM_cpu[i_m], nSpTcompPerSpM_Max);    
    nSpMcomp++;
  }

  CudaScalar<int>  nSpMcomp_cuda(&nSpMcomp);
  CudaScalar<int>  nSpBcompPerSpM_Max_cuda(&nSpBcompPerSpM_Max);
  CudaScalar<int>  nSpTcompPerSpM_Max_cuda(&nSpTcompPerSpM_Max);

  // Get the index of compatible hits
  auto mIndex_cpu = new int[nSpMcomp];  
  CpuMatrix< int > bIndex_cpu(nSpBcompPerSpM_Max, nSpMcomp); bIndex_cpu.Zeros();
  CpuMatrix< int > tIndex_cpu(nSpTcompPerSpM_Max, nSpMcomp); tIndex_cpu.Zeros();
    
  int mcnt(0);
  for (int i_m=0; i_m<nSpM; i_m++){
    if ( !nSpBcompPerSpM_cpu[i_m] || !nSpTcompPerSpM_cpu[i_m] ) continue;

    auto BottomCol = isCompatBottomMat_cpu.Get(0,i_m);
    int bcnt(0);
    for (int i=0; i<nSpB; i++){
      if (BottomCol[i]){
	bIndex_cpu.Set(bcnt,mcnt,i);
	bcnt++;
      }
    }
    
    auto TopCol = isCompatTopMat_cpu.Get(0,i_m);    
    int tcnt(0);
    for (int i=0; i<nSpT; i++){
      if (TopCol[i]){
	tIndex_cpu.Set(tcnt,mcnt,i);
	tcnt++;
      }
    }

    mIndex_cpu[mcnt] = i_m;    
    mcnt++;    
  }
  
  CudaVector< int > mIndex_cuda(nSpMcomp,    mIndex_cpu);
  CudaMatrix< int > bIndex_cuda(nSpBcompPerSpM_Max, nSpMcomp, &bIndex_cpu);
  CudaMatrix< int > tIndex_cuda(nSpTcompPerSpM_Max, nSpMcomp, &tIndex_cpu);

  /* -----------------------------------------
     Algorithm 1.B Matrix reduction
  -------------------------------------------*/
  
  CudaMatrix<float> spBcompMat_cuda(nSpBcompPerSpM_Max, nSpMcomp*6);
  CudaMatrix<float> spTcompMat_cuda(nSpTcompPerSpM_Max, nSpMcomp*6);

  dim3 RM_GridSize(nSpMcomp,1,1);
  dim3 RM_BlockSize; 

  // For bottom hits
  int offsetRM;
  offsetRM=0;
  while(offsetRM<nSpBcompPerSpM_Max){
    RM_BlockSize = dim3(fmin(m_config.maxBlockSize, nSpBcompPerSpM_Max-offsetRM),1,1);
    SeedfinderCudaKernels::reduceMatrix(RM_GridSize, RM_BlockSize,
					nSpB_cuda.Get(),
					spBmat_cuda.Get(0,0),
					nSpBcompPerSpM_Max_cuda.Get(),
					bIndex_cuda.Get(offsetRM,0),
					spBcompMat_cuda.Get(offsetRM,0)
					);
    offsetRM+=RM_BlockSize.x;    
  }

  // For top hits
  offsetRM=0;
  while(offsetRM<nSpTcompPerSpM_Max){
    RM_BlockSize = dim3(fmin(m_config.maxBlockSize, nSpTcompPerSpM_Max-offsetRM),1,1);
    SeedfinderCudaKernels::reduceMatrix(RM_GridSize, RM_BlockSize,
					nSpT_cuda.Get(),
					spTmat_cuda.Get(0,0),
					nSpTcompPerSpM_Max_cuda.Get(),
					tIndex_cuda.Get(offsetRM,0),
					spTcompMat_cuda.Get(offsetRM,0)
					);
    offsetRM+=RM_BlockSize.x;    
  }

  // For middle hits
  CpuMatrix<float>  spMcompMat_cpu (nSpMcomp, 6);
  CudaMatrix<float> spMcompMat_cuda(nSpMcomp, 6);
  for (int i_m=0; i_m<nSpMcomp; i_m++){
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

  CudaMatrix<float> circBcompMat_cuda(nSpBcompPerSpM_Max, nSpMcomp*6);  
  CudaMatrix<float> circTcompMat_cuda(nSpTcompPerSpM_Max, nSpMcomp*6);
    
  dim3 TC_GridSize(nSpMcomp,1,1);
  dim3 TC_BlockSize;
  
  // For bottom-middle doublets
  int offsetTC;
  offsetTC=0;
  while(offsetTC<nSpBcompPerSpM_Max){
    TC_BlockSize = dim3(fmin(m_config.maxBlockSize, nSpBcompPerSpM_Max-offsetTC),1,1);    
    SeedfinderCudaKernels::transformCoordinates(TC_GridSize,
						TC_BlockSize,
						true_cuda.Get(),
						spMcompMat_cuda.Get(0,0),
						nSpBcompPerSpM_Max_cuda.Get(),
						spBcompMat_cuda.Get(offsetTC,0),
						circBcompMat_cuda.Get(offsetTC,0));    
    offsetTC+=TC_BlockSize.x;
  }
  
  // For middle-top doublets
  offsetTC=0;
  while(offsetTC<nSpTcompPerSpM_Max){
    TC_BlockSize = dim3(fmin(m_config.maxBlockSize, nSpTcompPerSpM_Max-offsetTC),1,1);
    SeedfinderCudaKernels::transformCoordinates(TC_GridSize,
						TC_BlockSize,
						false_cuda.Get(),
						spMcompMat_cuda.Get(0,0),
						nSpTcompPerSpM_Max_cuda.Get(),
						spTcompMat_cuda.Get(offsetTC,0),
						circTcompMat_cuda.Get(offsetTC,0));    
    offsetTC+=TC_BlockSize.x;
  }
    
  // retreive middle-bottom doublet circ information
  CpuMatrix<float> circBcompMat_cpu(nSpBcompPerSpM_Max, nSpMcomp*6, &circBcompMat_cuda);
  
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

  CudaMatrix<int>   nTopPass_cuda(nSpBcompPerSpM_Max, nSpMcomp); nTopPass_cuda.Zeros();
  CudaMatrix<int>   tPassIndex_cuda(nTopPassLimit, nSpBcompPerSpM_Max*nSpMcomp); 
  CudaMatrix<float> curvatures_cuda(nTopPassLimit, nSpBcompPerSpM_Max*nSpMcomp);       
  CudaMatrix<float> impactparameters_cuda(nTopPassLimit, nSpBcompPerSpM_Max*nSpMcomp);

  CpuMatrix<int>    nTopPass_cpu(nSpBcompPerSpM_Max, nSpMcomp, true);
  CpuMatrix<int>    tPassIndex_cpu(nTopPassLimit, nSpBcompPerSpM_Max*nSpMcomp, true);
  CpuMatrix<float>  curvatures_cpu(nTopPassLimit, nSpBcompPerSpM_Max*nSpMcomp, true);
  CpuMatrix<float>  impactparameters_cpu(nTopPassLimit, nSpBcompPerSpM_Max*nSpMcomp, true);
  
  cudaStream_t cuStream;
  cudaStreamCreate(&cuStream);
   
  for (int i_m=0; i_m<=nSpMcomp; i_m++){    
    // -----------------------------------
    //  Algorithm 3. Triplet Search (TS)
    //------------------------------------

    cudaStreamSynchronize(cuStream);
    
    std::vector<std::pair<
      float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>> seedsPerSpM;
    
    std::vector<const InternalSpacePoint<external_spacepoint_t> *> tVec;
    std::vector<float> curvatures;
    std::vector<float> impactParameters;
    
    if (i_m < nSpMcomp){
          
      // For triplets collected at the previous iteration      
      int mIndex = mIndex_cpu[i_m];      
      dim3 TS_GridSize(nSpBcompPerSpM_cpu[mIndex],1,1);
      dim3 TS_BlockSize;
      
      int i_ts(0);    
      while ( offsetVec[i_ts] < nSpTcompPerSpM_cpu[mIndex] ){	
	TS_BlockSize = dim3(fmin(m_config.maxBlockSize, nSpTcompPerSpM_cpu[mIndex]-offsetVec[i_ts] ), 1,1); 
	SeedfinderCudaKernels::searchTriplet(TS_GridSize, TS_BlockSize,
					     offsetVec_cuda.Get(i_ts),
					     nSpMcomp_cuda.Get(),
					     spMcompMat_cuda.Get(i_m,0),
					     nSpBcompPerSpM_Max_cuda.Get(),
					     spBcompMat_cuda.Get(0,6*i_m),
					     nSpTcompPerSpM_Max_cuda.Get(),
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
					     tPassIndex_cuda.Get(0,nSpBcompPerSpM_Max*i_m),
					     curvatures_cuda.Get(0,nSpBcompPerSpM_Max*i_m),
					     impactparameters_cuda.Get(0,nSpBcompPerSpM_Max*i_m),
					     &cuStream
					     );
	i_ts++;
      }
      
      nTopPass_cpu.CopyD2H(nTopPass_cuda.Get(0,i_m),
			   nSpBcompPerSpM_Max,
			   nSpBcompPerSpM_Max*i_m,
			   &cuStream);
      tPassIndex_cpu.CopyD2H(tPassIndex_cuda.Get(0,nSpBcompPerSpM_Max*i_m),
			     nTopPassLimit*nSpBcompPerSpM_Max,
			     nTopPassLimit*nSpBcompPerSpM_Max*i_m,
			     &cuStream);
      
      curvatures_cpu.CopyD2H(curvatures_cuda.Get(0,nSpBcompPerSpM_Max*i_m),
			     nTopPassLimit*nSpBcompPerSpM_Max,
			     nTopPassLimit*nSpBcompPerSpM_Max*i_m,
			     &cuStream); 
      impactparameters_cpu.CopyD2H(impactparameters_cuda.Get(0,nSpBcompPerSpM_Max*i_m),
				   nTopPassLimit*nSpBcompPerSpM_Max,
				   nTopPassLimit*nSpBcompPerSpM_Max*i_m,
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
      
      for (int i_b=0; i_b<nSpBcompPerSpM_cpu[mIndex]; i_b++){
	int nTpass = *(nTopPass_cpu.Get(i_b,curID));
	if (nTpass==0) continue;	
	
	tVec.clear();
	curvatures.clear();
	impactParameters.clear();      
	float Zob = *(circBcompMat_cpu.Get(i_b,(curID)*6));

	std::vector< std::tuple< int, int, int > > indexVec;
	for(int i_t=0; i_t<nTpass; i_t++){
	  int g_tIndex = tIndexVec[*tPassIndex_cpu.Get(i_t,i_b+(curID)*nSpBcompPerSpM_Max)];
	  indexVec.push_back(std::make_tuple(g_tIndex,i_t,i_b));
	}
	sort(indexVec.begin(), indexVec.end()); 
	
	for(auto el: indexVec){
	  auto g_tIndex = std::get<0>(el);
	  auto tId      = std::get<1>(el);
	  auto bId      = std::get<2>(el);
	  
	  tVec.push_back(topSPvec[g_tIndex]);
	  curvatures.push_back(*curvatures_cpu.Get(tId,bId+(curID)*nSpBcompPerSpM_Max));
	  impactParameters.push_back(*impactparameters_cpu.Get(tId,bId+(curID)*nSpBcompPerSpM_Max));
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

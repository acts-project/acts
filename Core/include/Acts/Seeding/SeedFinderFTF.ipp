// SeedFinderFTF.ipp
// basing on ortho seed finder

#include "Acts/Definitions/Algebra.hpp"  //for M_PI
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SeedFinderFTFConfig.hpp"
#include "Acts/Seeding/SeedFinderUtils.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <type_traits>
#include <vector>

// core so in ACTS namespace

namespace Acts {

template <typename external_spacepoint_t>
SeedFinderFTF<external_spacepoint_t>::SeedFinderFTF(
    const SeedFinderFTFConfig<external_spacepoint_t> &config,
    const TrigFTF_GNN_Geometry<external_spacepoint_t> &GNNgeo)
    : m_config(config) {
  m_storage = new TrigFTF_GNN_DataStorage(GNNgeo);
}

template <typename external_spacepoint_t>
SeedFinderFTF<external_spacepoint_t>::~SeedFinderFTF() {
  delete m_storage;

  m_storage = nullptr;
}

// define loadspace points funciton
// when calling put input of vector<simspacepoints>, now can call space_point_t
template <typename external_spacepoint_t>
void SeedFinderFTF<external_spacepoint_t>::loadSpacePoints(
    const std::vector<FTF_SP<external_spacepoint_t>> &FTF_SP_vect) {
  for (auto &FTF_sp : FTF_SP_vect) {
    bool is_Pixel;
    if (FTF_sp.SP->sourceLinks().size() == 1) {  // pixels have 1 SL
      is_Pixel = true;
    } else {
      is_Pixel = false;
    }
    if (!is_Pixel)
      continue;

    m_storage->addSpacePoint(
        FTF_sp, (m_config.m_useClusterWidth > 0));  // add is a funciton FTFtype
  }

  m_config.m_phiSliceWidth = 2 * M_PI / m_config.m_nMaxPhiSlice;

  m_storage->sortByPhi();

  m_storage->generatePhiIndexing(1.5 * m_config.m_phiSliceWidth);
}


template <typename external_spacepoint_t>
void SeedFinderFTF<external_spacepoint_t>::runGNN_TrackFinder(std::vector<GNN_TrigTracklet<external_spacepoint_t>>&){

  // const int MaxEdges = 2000000;

  // const float cut_dphi_max      = 0.012;
  // const float cut_dcurv_max     = 0.001;
  // const float cut_tau_ratio_max = 0.007;
  // const float min_z0            = -2800; //roiDescriptor->zedMinus(); //blank for now, get from config eventually 
  // const float max_z0            = 2800; //roiDescriptor->zedPlus();

  // const float maxOuterRadius    = 550.0;
  // const float cut_zMinU = min_z0 // + maxOuterRadius*roiDescriptor->dzdrMinus();
  // const float cut_zMaxU = max_z0 // + maxOuterRadius*roiDescriptor->dzdrPlus();

  // float m_minR_squ = 1 ; //set earlier 
  // float m_maxCurv = 1 ; 

  // const float maxKappa_high_eta          = 0.8/m_minR_squ;
  // const float maxKappa_low_eta           = 0.6/m_minR_squ;

  // //1. loop over stages

  // int currentStage = 0;

  // const FASTRACK_CONNECTOR& conn = *(m_settings.m_conn); //turn into class name and template! 

  // std::vector<TrigFTF_GNN_Edge<external_spacepoint_t>> edgeStorage;
  
  // edgeStorage.reserve(MaxEdges);
  
  // int nEdges = 0;

}


template <typename external_spacepoint_t>
void SeedFinderFTF<external_spacepoint_t>::createSeeds(){
  
  std::vector<GNN_TrigTracklet<external_spacepoint_t>> vTracks; //make empty vector 

  vTracks.reserve(5000);

  runGNN_TrackFinder(vTracks); //returns filled vector 

  if(vTracks.empty()) return;

  m_triplets.clear(); //member of class , saying not declared, maybe public? 

  for(auto& track : vTracks) {
    for(auto& seed : track.m_seeds) { //access mmeber of GNN_TrigTracklet

      float newQ = seed.Q(); //function of TrigInDetTriplet
      if (m_config.m_LRTmode) { //dont have settings, just bool for now? 
	// In LRT mode penalize pixels in Triplets
	if(seed.s1().isPixel()) newQ+=1000; //functions of TrigSiSpacePointBase
	if(seed.s2().isPixel()) newQ+=1000;
	if(seed.s3().isPixel()) newQ+=1000;
      } else {
	// In normal (non LRT) mode penalise SSS by 1000, PSS (if enabled) and PPS by 10000
	if(seed.s3().isSCT()) { //functions of TrigSiSpacePointBase
	  newQ += seed.s1().isSCT() ? 1000.0 : 10000.0;
	} 
      }
      seed.Q(newQ);
      m_triplets.emplace_back(seed);
    }
  }
  vTracks.clear();


  
}


// // still to be developed
template <typename external_spacepoint_t>
template <typename input_container_t, typename output_container_t,
          typename callable_t>
void SeedFinderFTF<external_spacepoint_t>::createSeeds_old(
    const Acts::SeedFinderOptions &options,
    const input_container_t &spacePoints, output_container_t &out_cont,
    callable_t &&extract_coordinates) const {}

template <typename external_spacepoint_t>
template <typename input_container_t, typename callable_t>
std::vector<Seed<external_spacepoint_t>>
SeedFinderFTF<external_spacepoint_t>::createSeeds_old(
    const Acts::SeedFinderOptions &options,
    const input_container_t &spacePoints,
    callable_t &&extract_coordinates) const {
  std::vector<seed_t> r;
  createSeeds_old(options, spacePoints, r,
              std::forward<callable_t>(extract_coordinates));
  return r;
}

}  // namespace Acts

//SeedFinderFTF.ipp
//basing on ortho seed finder 

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SeedFinderFTFConfig.hpp"
#include "Acts/Seeding/SeedFinderUtils.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Definitions/Algebra.hpp" //for M_PI 

#include <fstream>
#include <vector>
#include <iostream>


#include <cmath>
#include <functional>
#include <numeric>
#include <type_traits>

//core so in ACTS namespace

namespace Acts {

template <typename external_spacepoint_t>
SeedFinderFTF<external_spacepoint_t>::SeedFinderFTF(
    const SeedFinderFTFConfig<external_spacepoint_t> &config, const TrigFTF_GNN_Geometry<external_spacepoint_t> &GNNgeo)
    : m_config(config) {
  

  // std::ifstream input_ifstream(m_config.fastrack_input_file.c_str(), std::ifstream::in) ;

  // FasTrackConnector input_fastrack(input_ifstream) ; 

  // TrigFTF_GNN_Geometry<external_spacepoint_t> mGNNgeo(m_config.input_vector, &input_fastrack);

  m_storage = new TrigFTF_GNN_DataStorage(GNNgeo);

 

}

template <typename external_spacepoint_t>
SeedFinderFTF<external_spacepoint_t>::~SeedFinderFTF(){

  delete m_storage ; 

  m_storage = nullptr ; 

}


//define loadspace points funciton 
 //when calling put input of vector<simspacepoints>, now can call space_point_t 
template <typename external_spacepoint_t>
void SeedFinderFTF<external_spacepoint_t>::loadSpacePoints(const std::vector<FTF_SP<external_spacepoint_t>>& FTF_SP_vect){
   
  for(auto& FTF_sp : FTF_SP_vect) {
    
    bool is_Pixel ; 
    if (FTF_sp.SP->sourceLinks().size() == 1) { //pixels have 1 SL
      is_Pixel = true ;
    }
    else {
      is_Pixel =  false ; 
    }
    // bool is_Pixel = FTF_sp.SP->isPixel(); //FTF actual object then sim is pointer 
    if(!is_Pixel) continue;

    m_storage->addSpacePoint(FTF_sp,(m_config.m_useClusterWidth > 0) ); //add is a funciton FTFtype 

  }

  m_config.m_phiSliceWidth = 2*M_PI/m_config.m_nMaxPhiSlice;

  m_storage->sortByPhi();

  m_storage->generatePhiIndexing(1.5*m_config.m_phiSliceWidth);

}


//still to be developed 
template <typename external_spacepoint_t>
template <typename input_container_t, typename output_container_t,
          typename callable_t>
void SeedFinderFTF<external_spacepoint_t>::createSeeds(
    const Acts::SeedFinderOptions &options,
    const input_container_t &spacePoints, output_container_t &out_cont,
    callable_t &&extract_coordinates) const {

} 

template <typename external_spacepoint_t>
template <typename input_container_t, typename callable_t>
std::vector<Seed<external_spacepoint_t>>
SeedFinderFTF<external_spacepoint_t>::createSeeds(
    const Acts::SeedFinderOptions &options,
    const input_container_t &spacePoints,
    callable_t &&extract_coordinates) const {
  std::vector<seed_t> r;
  createSeeds(options, spacePoints, r,
              std::forward<callable_t>(extract_coordinates));
  return r;
}



} //end of Acts namespace 

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
    const SeedFinderFTFConfig<external_spacepoint_t> &config)
    : m_config(config) {
  
  //checks if internal units funciton used 
  // if (not config.isInInternalUnits) {
  //   throw std::runtime_error(
  //       "SeedFinderOrthogonalConfig not in ACTS internal units in "
  //       "SeedFinderOrthogonal");

  // std::cout << "in FTF finder constructor 1 " << std::endl ;     

  std::ifstream input_ifstream(m_config.fastrack_input_file.c_str(), std::ifstream::in) ;
  // std::cout << "in FTF finder constructor 2 " << std::endl ;     

  FasTrackConnector input_fastrack(input_ifstream) ; 
  // std::cout << "in FTF finder constructor 3" << std::endl ;     

  TrigFTF_GNN_Geometry<external_spacepoint_t> mGNNgeo(m_config.input_vector, &input_fastrack);
  // std::cout << "in FTF finder constructor 4 " << std::endl ;     

  m_storage = new TrigFTF_GNN_DataStorage(mGNNgeo);
  // std::cout << "in FTF finder constructor 5 " << std::endl ;     

 

}

template <typename external_spacepoint_t>
SeedFinderFTF<external_spacepoint_t>::~SeedFinderFTF(){
  // std::cout << "in FTF finder destructor 1 " << std::endl ;     

  delete m_storage ; 
  // std::cout << "in FTF finder destructor 2 " << std::endl ;     

  m_storage = nullptr ; 
  // std::cout << "in FTF finder destructor 3 " << std::endl ;     

}


//define loadspace points funciton 
 //when calling put input of vector<simspacepoints>, now can call space_point_t 
template <typename external_spacepoint_t>
void SeedFinderFTF<external_spacepoint_t>::loadSpacePoints(const std::vector<FTF_SP<external_spacepoint_t>>& FTF_SP_vect){
  // std::cout << "start of lsp" << std::endl ;     
   
  for(auto& FTF_sp : FTF_SP_vect) {
    // std::cout << "in lsp 2 " << std::endl ;     

    bool is_Pixel = FTF_sp.SP->isPixel(); //FTF actual object then sim is pointer 
    if(!is_Pixel) continue;
    // std::cout << "in lsp 3 " << std::endl ;     

    //crashing here 
    m_storage->addSpacePoint(FTF_sp,(m_config.m_useTrigSeedML > 0) ); //add is a funciton FTFtype 
    // std::cout << "in lsp 4 " << std::endl ;     

  }
  // std::cout << "in lsp 5 " << std::endl ;     


  m_config.m_phiSliceWidth = 2*M_PI/m_config.m_nMaxPhiSlice;
  // std::cout << "in lsp 6 " << std::endl ;     

  m_storage->sortByPhi();
  // std::cout << "in lsp 7 " << std::endl ;     

  m_storage->generatePhiIndexing(1.5*m_config.m_phiSliceWidth);

  // std::cout << "end of lsp " << std::endl ;     

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

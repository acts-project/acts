//SeedFinderFTF.ipp
//accidently based on ipp, not sure if I'll neded 

//basing on ortho seed finder 

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SeedFinderFTFConfig.hpp"
#include "Acts/Seeding/SeedFinderUtils.hpp"
#include "Acts/Utilities/BinningType.hpp"

//for load space points 
// #include "Acts/Seeding/FastTrackConnector.h"
#include "Acts/Seeding/GNN_DataStorage.h"
// #include "Acts/Seeding/SeedFinderFTF.hpp"


#include <cmath>
#include <functional>
#include <numeric>
#include <type_traits>

//core so in ACTS namespace

namespace Acts {

//valid tuple ortho range function RH and LH 
//valid tuple funtion 


//this is the definition of constructos 

template <typename external_spacepoint_t>
SeedFinderFTF<external_spacepoint_t>::SeedFinderFTF(
    const SeedFinderFTFConfig<external_spacepoint_t> &config)
    : m_config(config) {
  m_storage = new TrigFTF_GNN_DataStorage<external_spacepoint_t>(m_config.mGNNgeo); //when define find relevant input 
  
  //schecks if internal units funciton used 
  // if (not config.isInInternalUnits) {
  //   throw std::runtime_error(
  //       "SeedFinderOrthogonalConfig not in ACTS internal units in "
  //       "SeedFinderOrthogonal");
}
//destructor too? deletes m_storage 


//define loadspace points funciton 
 //when calling put input of vector<simspacepoints>, now can call space_point_t 
template <typename external_spacepoint_t>
void SeedFinderFTF<external_spacepoint_t>::loadSpacePoints(const std::vector<external_spacepoint_t>& SP){ 

  for(auto sp : SP) {
    //could check if pixel as pixels only have 1 source link (strip have 2)
    bool isPixel = sp.isPixel(); 

    if(!isPixel) continue;
    // //think not using trigseedML for now 
    // //when called input should be simspace point 
    
    m_storage->addSpacePoint((sp), (m_config.m_useTrigSeedML > 0));

    std::cout<<("in seed finder load space point function") ; 
  }

  //cant find M_PI set to number for now 
  int M_PI = 1 ; 
  m_phiSliceWidth = 2*M_PI/m_config.m_nMaxPhiSlice;


  m_storage->sortByPhi();
  m_storage->generatePhiIndexing(1.5*m_phiSliceWidth);

}
//filter canditates function 
//fucntion processFromMiddleSP 

//create tree function- called in seeds 

//need create seeds function- at same level as TrigTrackSeedGenerator_itk  
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
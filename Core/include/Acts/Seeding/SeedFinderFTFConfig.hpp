//basing on ortho config
#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Seeding/GNN_DataStorage.hpp"

#include <memory> 

//core algorithm so in acts namespace 
namespace Acts { 

template <typename T>
class SeedFilter;


template <typename SpacePoint> 
struct SeedFinderFTFConfig {

  // // how many sigmas of scattering angle should be considered?
  float sigmaScattering = 5;

  //Seed cut 
  float minPt = 400. * Acts::UnitConstants::MeV;

  ///////////some declared not filled in by reco: //////
  std::shared_ptr<Acts::SeedFilter<SpacePoint>> seedFilter;

  //   //detector ROI   
  //   // derived values, set on SeedFinder construction
  float highland = 0; 
  float maxScatteringAngle2 = 0;
  //   bool isInInternalUnits = false;
  ///for load space points 
  unsigned int maxSeedsPerSpM = 5;

  float m_phiSliceWidth ; 
  float m_nMaxPhiSlice ; 
  bool m_useClusterWidth = false ; 
  std::string fastrack_input_file ; 
  std::vector<TrigInDetSiLayer> input_vector ;

  ////
  //2 member functions 
  SeedFinderFTFConfig calculateDerivedQuantities() const {
    //thorw statement if the isInternalUnits member is false, ie if dont call this function 
    SeedFinderFTFConfig config = *this;
    //use a formula to calculate scattering 

    return config;
  }

  SeedFinderFTFConfig toInternalUnits() const {
   //throw statement if the isInternalUnits member is false, ie if dont call this function 
   SeedFinderFTFConfig config = *this;
   //devides inputs by 1mm, all ones input  
   //changes memeber inInInternalUnits to true 
    return config;
  }


}; //end of config struct 

} //end of namespace ACTS 
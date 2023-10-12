// basing on ortho config
#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
// #include "Acts/Seeding/GNN_DataStorage.hpp"
#include "Acts/Seeding/TrigBase.hpp" //definition of Trigsispacepoint base and trigtriplets 



#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"

#include <memory>

// core algorithm so in acts namespace
namespace Acts {

template <typename T> class SeedFilter;

template <typename SpacePoint> struct SeedFinderFTFConfig {

  // // how many sigmas of scattering angle should be considered?
  float sigmaScattering = 5;

  // Seed cut
  float minPt = 400. * Acts::UnitConstants::MeV;

  ///////////some declared not filled in by reco: //////
  std::shared_ptr<Acts::SeedFilter<SpacePoint>> seedFilter;

  //   //detector ROI
  //   // derived values, set on SeedFinder construction
  float highland = 0;
  float maxScatteringAngle2 = 0;
  //   bool isInInternalUnits = false;
  /// for load space points
  unsigned int maxSeedsPerSpM = 5;

  float m_phiSliceWidth;
  float m_nMaxPhiSlice;
  bool m_useClusterWidth = false;
  std::string fastrack_input_file;
  std::vector<TrigInDetSiLayer> m_layerGeometry; //used to be input_vector 

  //for run funciton 
  //m_settings: 
  bool m_LRTmode = true ; //eventually want to set from full chaing 
  bool m_useEtaBinning = true ; 
  bool m_doubletFilterRZ = true ; 
  float m_minDeltaRadius = 5.0 ; //eventaully set in config or to equivelant acts 2.0 but increasing to test loops 
  // float m_maxDeltaRadius = 270.0 ; 
  float m_tripletD0Max = 4.0 ; //m_settings 
  unsigned int m_maxTripletBufferLength = 3 ;  

  //ROI: 
  bool containsPhi() {
    return false; 
    //need to implement this funciton 
  }

  ////
  // 2 member functions
  SeedFinderFTFConfig calculateDerivedQuantities() const {
    // thorw statement if the isInternalUnits member is false, ie if dont call
    // this function
    SeedFinderFTFConfig config = *this;
    // use a formula to calculate scattering

    return config;
  }

  SeedFinderFTFConfig toInternalUnits() const {
    // throw statement if the isInternalUnits member is false, ie if dont call
    // this function
    SeedFinderFTFConfig config = *this;
    // devides inputs by 1mm, all ones input
    // changes memeber inInInternalUnits to true
    return config;
  }

}; // end of config struct

} // namespace Acts
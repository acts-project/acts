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



  // limiting location of measurements  
  float rMin = 33 * Acts::UnitConstants::mm;
  float rMax = 600 * Acts::UnitConstants::mm; 

  // Seed Cuts
  // minimum distance in r between middle and top SP in one seed
  float deltaRMinTopSP = 5 * Acts::UnitConstants::mm;
  // maximum distance in r between middle and top SP in one seed
  float deltaRMaxTopSP = 270 * Acts::UnitConstants::mm;
  // minimum distance in r between middle and bottom SP in one seed
  float deltaRMinBottomSP = 5 * Acts::UnitConstants::mm;
  // maximum distance in r between middle and bottom SP in one seed
  float deltaRMaxBottomSP = 270 * Acts::UnitConstants::mm;

  //Geometry, detector ROI, limitng collisioon regions
  float collisionRegionMin = -150 * Acts::UnitConstants::mm;
  float collisionRegionMax = +150 * Acts::UnitConstants::mm;
  float zMin = -2800 * Acts::UnitConstants::mm;
  float zMax = 2800 * Acts::UnitConstants::mm;

  // for how many seeds can one SpacePoint be the middle SpacePoint?
  unsigned int maxSeedsPerSpM = 5;

  //Seed cuts 
  // cot of maximum theta angle, equivalent to 2.7 eta (pseudorapidity)
  float cotThetaMax = 7.40627;

  // how many sigmas of scattering angle should be considered?
  float sigmaScattering = 5;

  // average radiation lengths of material on the length of a seed. used for
  float radLengthPerSeed = 0.05; 

  //Seed cut 
  float minPt = 400. * Acts::UnitConstants::MeV;

  // impact parameter
  float impactMax = 20. * Acts::UnitConstants::mm;

  // enable cut on the compatibility between interaction point and SPs
  bool interactionPointCut = false;

  // cut to the maximum value of delta z between SPs
  float deltaZMax = std::numeric_limits<float>::infinity() * Acts::UnitConstants::mm;

  // Upper pt limit for scattering calculation
  float maxPtScattering = 10 * Acts::UnitConstants::GeV;

  // radial range for middle SP
  // range defined in vector for each z region
  std::vector<std::vector<float>> rRangeMiddleSP; 
  // variable range based on SP radius
  bool useVariableMiddleSPRange = true;

  // seed confirmation
  bool seedConfirmation = false;
  // parameters for central seed confirmation
  SeedConfirmationRangeConfig centralSeedConfirmationRange;
  // parameters for forward seed confirmation
  SeedConfirmationRangeConfig forwardSeedConfirmationRange;

  ///////////some declared not filled in by reco: //////
  std::shared_ptr<Acts::SeedFilter<SpacePoint>> seedFilter;

  //   //detector ROI   
  float phiMin = -M_PI;
  float phiMax = M_PI;
  //    // radial range for middle SP
  //   float deltaRMiddleMinSPRange = 10. * Acts::UnitConstants::mm;
  //   float deltaRMiddleMaxSPRange = 10. * Acts::UnitConstants::mm;
  //   float rMinMiddle = 60.f * Acts::UnitConstants::mm;
  //   float rMaxMiddle = 120.f * Acts::UnitConstants::mm;
  //   float deltaPhiMax = 0.085; 
  //   // skip top SPs based on cotTheta sorting when producing triplets
  //   bool skipPreviousTopSP = false;
  //   // derived values, set on SeedFinder construction
  float highland = 0; 
  float maxScatteringAngle2 = 0;
  //   bool isInInternalUnits = false;


  ///for load space points 

  float m_phiSliceWidth ; 
  float m_nMaxPhiSlice ; 
  bool m_useTrigSeedML = false ; 
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
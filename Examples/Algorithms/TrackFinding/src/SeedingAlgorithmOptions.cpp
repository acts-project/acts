// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/SeedingAlgorithmOptions.hpp"

#include <string>

#include <boost/program_options.hpp>

void ActsExamples::Options::addSeedingAlgorithmOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();

  opt("trkseedfilter-deltaInvHelixDiameter", value<double>()->default_value(0.00003),
      "the allowed delta between two inverted seed radii for them to be considered compatible. delta of inverse radius leads to this value being the cutoff. unit is 1/mm. default value of 0.00003 leads to all helices with radius>33m to be considered compatible (in mm)");

  opt("trkseedfilter-compatSeedWeight", value<double>()->default_value(200),
      "seed weight increased by this value if a compatible seed has been found");

  opt("trkseedfilter-compatSeedLimit", value<size_t>()->default_value(2),
      "number of compatible seeds considered for weight increase");

  opt("trkseedfinder-minPt", value<double>()->default_value(400.),
      "lower pt cutoff for seeds (in MeV)");

  opt("trkseedfinder-deltaRMin", value<double>()->default_value(5),
      "minimum distance in r between two measurements within one seed (in mm)");

  opt("trkseedfinder-deltaRMax", value<double>()->default_value(270),
      "maximum distance in r between two measurements within one seed (in mm)");

  opt("trkseedfinder-impactMax", value<double>()->default_value(20),
      "max impact parameter (in mm)");

  opt("trkseedfinder-sigmaScattering", value<double>()->default_value(5),
      "number of sigmas of scattering angle considered");

  opt("trkseedfinder-maxPtScattering", value<double>()->default_value(10),
      "Upper pt limit for scattering calculation (in GeV)");

  opt("trkseedfinder-maxSeedsPerSpM ", value<unsigned int>()->default_value(5),
      "for how many seeds can one SpacePoint be the middle SpacePoint");

  opt("trkseedfinder-collisionRegionMin", value<double>()->default_value(-150),
      "min collision region in z (in mm)");

  opt("trkseedfinder-collisionRegionMax", value<double>()->default_value(150),
      "max collision region in z (in mm)");

  opt("trkseedfinder-zMin", value<double>()->default_value(-2800),
      "min z (in mm)");

  opt("trkseedfinder-zMax", value<double>()->default_value(2800),
      "max z (in mm)");

  opt("trkseedfinder-rMax", value<double>()->default_value(600),
      "max r (in mm)");

  opt("trkseedfinder-rMin", value<double>()->default_value(33),
      "min r, if rMin is smaller than impactMax, the bin size will be 2*pi, which will make seeding very slow! (in mm)");

  opt("trkseedfinder-radLengthPerSeed", value<double>()->default_value(0.05),
      "average radiation lengths of material on the length of a seed. used for scattering.");

}

void ActsExamples::Options::addMLOutputOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();

  opt("output-ML", value<bool>()->default_value(false), "true if output should be ML friendly for EA algorithm");

}

ActsExamples::SeedingAlgorithm::Config
ActsExamples::Options::readSeedingAlgorithmConfig(
    const ActsExamples::Options::Variables& vars) {
  using namespace ActsExamples;
  using namespace Acts::UnitConstants;

  SeedingAlgorithm::Config cfg;

  cfg.seedFilterConfig.deltaInvHelixDiameter = vars["trkseedfilter-deltaInvHelixDiameter"].as<double>() * 1. / Acts::UnitConstants::mm; 
  cfg.seedFilterConfig.compatSeedWeight = vars["trkseedfilter-compatSeedWeight"].as<double>();
  cfg.seedFilterConfig.deltaRMin = vars["trkseedfinder-deltaRMin"].as<double>() * Acts::UnitConstants::mm;
  cfg.seedFilterConfig.maxSeedsPerSpM = vars["trkseedfinder-maxSeedsPerSpM"].as<unsigned int>();
  cfg.seedFilterConfig.compatSeedLimit = vars["trkseedfilter-compatSeedLimit"].as<size_t>();   

  cfg.seedFinderConfig.minPt = vars["trkseedfinder-minPt"].as<double>() * Acts::UnitConstants::MeV;
  cfg.seedFinderConfig.deltaRMin = vars["trkseedfinder-deltaRMin"].as<double>() * Acts::UnitConstants::mm;
  cfg.seedFinderConfig.deltaRMax = vars["trkseedfinder-deltaRMax"].as<double>() * Acts::UnitConstants::mm;
  cfg.seedFinderConfig.impactMax = vars["trkseedfinder-impactMax"].as<double>() * Acts::UnitConstants::mm;
  cfg.seedFinderConfig.sigmaScattering = vars["trkseedfinder-sigmaScattering"].as<double>();
  cfg.seedFinderConfig.maxPtScattering = vars["trkseedfinder-maxPtScattering"].as<double>() * Acts::UnitConstants::GeV;
  cfg.seedFinderConfig.maxSeedsPerSpM = vars["trkseedfinder-maxSeedsPerSpM"].as<unsigned int>();
  cfg.seedFinderConfig.collisionRegionMin = vars["trkseedfinder-collisionRegionMin"].as<double>() * Acts::UnitConstants::mm;
  cfg.seedFinderConfig.collisionRegionMax = vars["trkseedfinder-collisionRegionMax"].as<double>() * Acts::UnitConstants::mm;
  cfg.seedFinderConfig.zMin = vars["trkseedfinder-zMin"].as<double>() * Acts::UnitConstants::mm;
  cfg.seedFinderConfig.zMax = vars["trkseedfinder-zMax"].as<double>() * Acts::UnitConstants::mm;
  cfg.seedFinderConfig.rMax = vars["trkseedfinder-rMax"].as<double>() * Acts::UnitConstants::mm;
  cfg.seedFinderConfig.rMin = vars["trkseedfinder-rMin"].as<double>() * Acts::UnitConstants::mm;
  cfg.seedFinderConfig.radLengthPerSeed = vars["trkseedfinder-radLengthPerSeed"].as<double>(); 

  cfg.gridConfig.minPt = vars["trkseedfinder-minPt"].as<double>() * Acts::UnitConstants::MeV;
  cfg.gridConfig.rMax = vars["trkseedfinder-rMax"].as<double>() * Acts::UnitConstants::mm;
  cfg.gridConfig.zMax = vars["trkseedfinder-zMax"].as<double>() * Acts::UnitConstants::mm;
  cfg.gridConfig.zMin = vars["trkseedfinder-zMin"].as<double>() * Acts::UnitConstants::mm;
  cfg.gridConfig.deltaRMax = vars["trkseedfinder-deltaRMax"].as<double>() * Acts::UnitConstants::mm;
  cfg.gridConfig.impactMax = vars["trkseedfinder-impactMax"].as<double>() * Acts::UnitConstants::mm;

  return cfg;
}

bool ActsExamples::Options::readMLOutputConfig(const ActsExamples::Options::Variables& vm) {
  bool outputIsML = false;
  if (vm.count("output-ML")) {
    outputIsML = vm["output-ML"].as<bool>();
  }
  return outputIsML;
}

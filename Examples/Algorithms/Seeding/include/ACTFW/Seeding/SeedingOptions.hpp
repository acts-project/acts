// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>

#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/Seeding/SeedingAlgorithm.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

namespace po = boost::program_options;
using namespace Acts::UnitLiterals;

namespace FW {

namespace Options {

/// @brief Seeding options : base
///
/// @tparam aopt_t Type of the options class from boost
template <typename aopt_t>
void addSeedingOptions(aopt_t& opt) {
  opt.add_options()("seed-collection",
                    po::value<std::string>()->default_value("seed-collection"),
                    "Name of the output seed collection.")(
      "seed-pt-range", po::value<Interval>()->value_name("MIN:MAX"),
      "Transverse momentum range [MeV] for the seed finding.")(
      "seed-eta-range", po::value<Interval>()->value_name("MIN:MAX"),
      "Pseudorapidity range for the seed finding.")(
      "seed-delta-r-range", po::value<Interval>()->value_name("MIN:MAX"),
      "Min/max distance [mm] between two measurements within one seed.")(
      "seed-d0-range", po::value<Interval>()->value_name("MIN:MAX"),
      "Min/max transverse impact parameter [mm] allowed for a seed.")(
      "seed-z0-range", po::value<Interval>()->value_name("MIN:MAX"),
      "Min/max longitudinal impact parameter [mm] allowed for a seed.")(
      "seed-sp-r-range", po::value<Interval>()->value_name("MIN:MAX"),
      "Min/max radius [mm] for space points considered for seeding.")(
      "seed-sp-z-range", po::value<Interval>()->value_name("MIN:MAX"),
      "Min/max z [mm] for space points considered for seeding.")(
      "seed-sp-phi-range", po::value<Interval>()->value_name("MIN:MAX"),
      "Min/max phi for space points considered for seeding.")(
      "seed-scattering-tX0", po::value<double>()->default_value(0.05),
      "Averge t/X0 seed by a seed for scattering inclusion")(
      "seed-scattering-nsigma", po::value<double>()->default_value(5.),
      "Number of sigmas of scattering component considered for a seed.")(
      "seed-alignment-r-nsigma", po::value<double>()->default_value(0.),
      "Number of sigmas of radial misalignment considered for a seed.")(
      "seed-alignment-z-nsigma", po::value<double>()->default_value(0.),
      "Number of sigmas of longitudinal mislaignment considered for a seed.")(
      "seed-nseeds-per-sp", po::value<int>()->default_value(5),
      "Allowed seed candidates per space point.");
}

/// @brief Seeding options : base
///
/// @tparam aopt_t Type of the options class from boost
template <typename aopt_t>
void addCudaSeedingOptions(aopt_t& opt) {
  opt.add_options()("seed-cuda-max-block-size",
                    po::value<int>()->default_value(1024),
                    "For CUDA implementation only: maximum block side.")(
      "seed-cuda-nseeds-per-sp-limit", po::value<int>()->default_value(100),
      "For CUDA implementation only: limit of number of seeds per spacepoint.")(
      "seed-cuda-nseeds-per-sp-avg", po::value<int>()->default_value(2),
      "For CUDA implementation only: average number of seeds per spacepoint.");
}

/// read the evgen options and return a Config file
///
/// @tparam vmap_t is the Type of the Parameter map to be read out
/// @tparam bfield_t is the Type of the Magnetic field
///
/// @param vm is the parameter map for the options
/// @param magField is the magnetic field objects as shared pointer
/// @param tGeometry is the tracking geometry object
///
/// @returns a Config object for the ExtrapolationAlgorithm
template <typename vmap_t>
typename FW::SeedingAlgorithm<Acts::SeedFinder<SimHit>>::Config
readSeedingConfig(const vmap_t& vm) {
  using SeedFinderConfig = Acts::SeedFinderConfig<SimHit>;
  using SeedFinder = Acts::SeedFinder<SimHit>;

  SeedFinderConfig seedFinderCfg;
  // TODO connect with options
  seedFinderCfg.minPt = 400_MeV;
  seedFinderCfg.cotThetaMax = 7.40627;
  seedFinderCfg.deltaRMin = 5_mm;
  seedFinderCfg.deltaRMax = 270_mm;
  seedFinderCfg.impactMax = 20_mm;
  seedFinderCfg.sigmaScattering = 5;
  seedFinderCfg.maxSeedsPerSpM = 5;
  seedFinderCfg.collisionRegionMin = -150_mm;
  seedFinderCfg.collisionRegionMax = +150_mm;
  seedFinderCfg.phiMin = -M_PI;
  seedFinderCfg.phiMax = M_PI;
  seedFinderCfg.zMin = -2800_mm;
  seedFinderCfg.zMax = 2800_mm;
  seedFinderCfg.rMax = 600_mm;
  seedFinderCfg.rMin = 33_mm;
  seedFinderCfg.bFieldInZ = 0.00208;
  seedFinderCfg.radLengthPerSeed = 0.05;
  seedFinderCfg.zAlign = 0;
  seedFinderCfg.rAlign = 0;
  seedFinderCfg.sigmaError = 5;

  SeedFinder seedFinder(seedFinderCfg);

  SeedingAlgorithm<SeedFinder>::Config seedAlgConfig(std::move(seedFinder));

  return seedAlgConfig;
  //
}

}  // namespace Options
}  // namespace FW

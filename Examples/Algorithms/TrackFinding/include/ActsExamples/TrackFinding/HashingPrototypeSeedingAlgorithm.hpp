// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding2/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding2/TripletSeeder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>

namespace ActsExamples {

/// Construct track seeds from space points.
class HashingPrototypeSeedingAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input space point collections.
    std::string inputSpacePoints;
    /// Output track seed collection.
    std::string outputSeeds;
    /// Output buckets collection.
    std::string outputBuckets;

    // General seeding parameters

    /// Magnetic field in z direction
    float bFieldInZ = 2 * Acts::UnitConstants::T;
    /// minimum pT
    float minPt = 0.4 * Acts::UnitConstants::GeV;
    /// maximum forward direction expressed as cot(theta)
    float cotThetaMax = 10.01788;  // equivalent to eta = 3 (pseudorapidity)
    /// maximum impact parameter in mm
    float impactMax = 20 * Acts::UnitConstants::mm;
    /// Minimum radial distance between two doublet components (prefer
    /// deltaRMinTop and deltaRMinBottom to set separate values for top and
    /// bottom space points)
    float deltaRMin = 5 * Acts::UnitConstants::mm;
    /// Maximum radial distance between two doublet components (prefer
    /// deltaRMaxTop and deltaRMacBottom to set separate values for top and
    /// bottom space points)
    float deltaRMax = 270 * Acts::UnitConstants::mm;
    /// Minimum radial distance between middle-top doublet components
    float deltaRMinTop = std::numeric_limits<float>::quiet_NaN();
    /// Maximum radial distance between middle-outer doublet components
    float deltaRMaxTop = std::numeric_limits<float>::quiet_NaN();
    /// Minimum radial distance between bottom-middle doublet components
    float deltaRMinBottom = std::numeric_limits<float>::quiet_NaN();
    /// Maximum radial distance between bottom-middle doublet components
    float deltaRMaxBottom = std::numeric_limits<float>::quiet_NaN();

    // Seeding parameters used to define the cuts on space-point doublets

    /// Minimal value of z-distance between space-points in doublet
    float deltaZMin = -std::numeric_limits<float>::infinity();
    /// Maximum value of z-distance between space-points in doublet
    float deltaZMax = std::numeric_limits<float>::infinity();

    // Seed finder doublet cuts

    /// Enable cut on the compatibility between interaction point and doublet,
    /// this is an useful approximation to speed up the seeding
    bool interactionPointCut = false;

    /// Limiting location of collision region in z-axis used to check if doublet
    /// origin is within reasonable bounds
    float collisionRegionMin = -150 * Acts::UnitConstants::mm;
    float collisionRegionMax = +150 * Acts::UnitConstants::mm;

    /// Parameter which can loosen the tolerance of the track seed to form a
    /// helix. This is useful for e.g. misaligned seeding.
    float helixCutTolerance = 1;

    // Seed finder triplet cuts

    /// Number of sigmas of scattering angle to be considered in the minimum pT
    /// scattering term
    float sigmaScattering = 5;
    /// Term that accounts for the thickness of scattering medium in radiation
    /// lengths in the Lynch & Dahl correction to the Highland equation default
    /// is 5%
    /// TODO: necessary to make amount of material dependent on detector region?
    float radLengthPerSeed = 0.05;

    /// Tolerance parameter used to check the compatibility of space-point
    /// coordinates in xyz. This is only used in a detector specific check for
    /// strip modules
    float toleranceParam = 1.1 * Acts::UnitConstants::mm;

    // Seed filter parameters

    /// Allowed difference in curvature (inverted seed radii) between two
    /// compatible seeds
    float deltaInvHelixDiameter = 0.00003 * (1 / Acts::UnitConstants::mm);
    /// Seed weight/score is increased by this value if a compatible seed has
    /// been found. This is the c1 factor in the seed score calculation (w = c1
    /// * Nt - c2 * d0 - c3 * z0)
    float compatSeedWeight = 200;
    /// The transverse impact parameters (d0) is multiplied by this factor and
    /// subtracted from weight. This is the c2 factor in the seed score
    /// calculation (w = c1 * Nt - c2 * d0 - c3 * z0)
    float impactWeightFactor = 1;
    /// The logitudinal impact parameters (z0) is multiplied by this factor and
    /// subtracted from weight. This is the c3 factor in the seed score
    /// calculation (w = c1 * Nt - c2 * d0 - c3 * z0)
    float zOriginWeightFactor = 1;
    /// Maximum number (minus one) of accepted seeds per middle space-point
    /// In dense environments many seeds may be found per middle space-point
    /// Only seeds with the highest weight will be kept if this limit is reached
    unsigned int maxSeedsPerSpM = 5;
    /// Maximum limit to number of compatible space-point used in score
    /// calculation. We increase by c1 the weight calculation for each
    /// compatible space-point until we reach compatSeedLimit
    std::size_t compatSeedLimit = 2;

    /// Increment in seed weight if the number of compatible seeds is larger
    /// than numSeedIncrement, this is used in case of high occupancy scenarios
    /// if we want to increase the weight of the seed by seedWeightIncrement
    /// when the number of compatible seeds is higher than a certain value
    float seedWeightIncrement = 0;
    float numSeedIncrement = std::numeric_limits<float>::infinity();

    /// Enable quality seed confirmation, this is different than default seeding
    /// confirmation because it can also be defined for different (r, z) regions
    /// of the detector (e.g. forward or central region) by
    /// SeedConfirmationRange. Seeds are classified as "high-quality" seeds and
    /// normal quality seeds. Normal quality seeds are only selected if no other
    /// "high-quality" seeds has been found for that inner-middle doublet.
    bool seedConfirmation = false;
    /// Contains parameters for central seed confirmation
    Acts::SeedConfirmationRangeConfig centralSeedConfirmationRange;
    /// Contains parameters for forward seed confirmation
    Acts::SeedConfirmationRangeConfig forwardSeedConfirmationRange;

    /// If seedConfirmation is true we classify seeds as "high-quality" seeds.
    /// Seeds that are not confirmed as "high-quality" are only selected if no
    /// other "high-quality" seed has been found for that inner-middle doublet
    /// Maximum number of normal seeds (not classified as "high-quality" seeds)
    /// in seed confirmation
    std::uint32_t maxSeedsPerSpMConf = 5;
    /// Maximum number of "high-quality" seeds for each inner-middle SP-dublet
    /// in seed confirmation. If the limit is reached we check if there is a
    /// normal quality seed to be replaced
    std::uint32_t maxQualitySeedsPerSpMConf = 5;

    /// Use deltaR between top and middle SP instead of top radius to search for
    /// compatible SPs
    bool useDeltaRinsteadOfTopRadius = false;

    // other

    /// Connect custom selections on the space points or to the doublet
    /// compatibility
    bool useExtraCuts = false;

    // hashing training

    /// Random seed for Annoy
    std::uint32_t annoySeed = 123456789;

    /// Number of features to use
    std::int32_t f = 1;

    // hashing inference

    /// Size of the buckets = number of space points in the bucket
    std::uint32_t bucketSize = 10;
    /// Number of zBins
    std::uint32_t zBins = 0;
    /// Number of phiBins
    std::uint32_t phiBins = 50;

    /// Layer selection
    double layerRMin = 25;
    double layerRMax = 40;
    double layerZMin = -550;
    double layerZMax = 550;
  };

  /// Construct the seeding algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  HashingPrototypeSeedingAlgorithm(
      Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Run the seeding algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  Acts::BroadTripletSeedFilter::Config m_filterConfig;
  std::unique_ptr<const Acts::Logger> m_filterLogger;
  std::optional<Acts::TripletSeeder> m_seedFinder;

  Acts::Delegate<bool(const SimSpacePoint&)> m_spacePointSelector{
      Acts::DelegateFuncTag<voidSpacePointSelector>{}};

  static bool voidSpacePointSelector(const SimSpacePoint& /*sp*/) {
    return true;
  }

  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this,
                                                            "InputSpacePoints"};

  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};
  WriteDataHandle<std::vector<SimSpacePointContainer>> m_outputBuckets{
      this, "OutputBuckets"};
};

}  // namespace ActsExamples

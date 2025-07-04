// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Seeding2/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding2/BroadTripletSeedFinder.hpp"
#include "Acts/Seeding2/CylindricalSpacePointGrid2.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace ActsExamples {

/// Construct track seeds from space points.
class GridTripletSeedingAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input space point collections.
    std::string inputSpacePoints;
    /// Output track seed collection.
    std::string outputSeeds;

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

    // Seeding parameters used in the space-point grid creation and bin finding

    /// minimum extension of sensitive detector layer relevant for seeding as
    /// distance from x=y=0 (i.e. in r)
    /// WARNING: if rMin is smaller than impactMax, the bin size will be 2*pi,
    /// which will make seeding very slow!
    float rMin = 0 * Acts::UnitConstants::mm;
    /// maximum extension of sensitive detector layer relevant for seeding as
    /// distance from x=y=0 (i.e. in r)
    float rMax = 600 * Acts::UnitConstants::mm;
    /// minimum extension of sensitive detector layer relevant for seeding in
    /// negative direction in z
    float zMin = -2800 * Acts::UnitConstants::mm;
    /// maximum extension of sensitive detector layer relevant for seeding in
    /// positive direction in z
    float zMax = 2800 * Acts::UnitConstants::mm;
    /// minimum phi value for phiAxis construction
    float phiMin = -std::numbers::pi_v<float>;
    /// maximum phi value for phiAxis construction
    float phiMax = std::numbers::pi_v<float>;
    /// Multiplicator for the number of phi-bins. The minimum number of phi-bins
    /// depends on min_pt, magnetic field: 2*pi/(minPT particle phi-deflection).
    /// phiBinDeflectionCoverage is a multiplier for this number. If
    /// numPhiNeighbors (in the configuration of the BinFinders) is configured
    /// to return 1 neighbor on either side of the current phi-bin (and you want
    /// to cover the full phi-range of minPT), leave this at 1.
    int phiBinDeflectionCoverage = 1;
    /// maximum number of phi bins
    int maxPhiBins = 10000;

    /// vector containing the map of z bins in the top and bottom layers
    std::vector<std::pair<int, int>> zBinNeighborsTop;
    std::vector<std::pair<int, int>> zBinNeighborsBottom;
    /// number of phiBin neighbors at each side of the current bin that will be
    /// used to search for SPs
    int numPhiNeighbors = 1;

    /// Vector containing the z-bin edges for non equidistant binning in z
    std::vector<float> zBinEdges;

    /// Order of z bins to loop over when searching for SPs
    std::vector<std::size_t> zBinsCustomLooping;

    // Seeding parameters used to define the region of interest for middle
    // space-point

    /// Radial range for middle space-point
    /// The range can be defined manually with (rMinMiddle, rMaxMiddle). If
    /// useVariableMiddleSPRange is set to false and the vector rRangeMiddleSP
    /// is empty, we use (rMinMiddle, rMaxMiddle) to cut the middle space-points
    float rMinMiddle = 60 * Acts::UnitConstants::mm;
    float rMaxMiddle = 120 * Acts::UnitConstants::mm;
    /// If useVariableMiddleSPRange is set to false, the vector rRangeMiddleSP
    /// can be used to define a fixed r range for each z bin: {{rMin, rMax},
    /// ...}
    bool useVariableMiddleSPRange = false;
    /// Range defined in vector for each z bin
    std::vector<std::vector<float>> rRangeMiddleSP;

    float deltaRMiddleMinSPRange = 10 * Acts::UnitConstants::mm;
    float deltaRMiddleMaxSPRange = 10 * Acts::UnitConstants::mm;

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
    /// Maximum transverse momentum for scattering calculation
    float maxPtScattering = 10 * Acts::UnitConstants::GeV;

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

    /// Use deltaR between top and middle SP instead of top radius to search for
    /// compatible SPs
    bool useDeltaRinsteadOfTopRadius = false;

    // Seeding parameters used for quality seed confirmation

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

    // other

    /// Connect custom selections on the space points or to the doublet
    /// compatibility
    bool useExtraCuts = false;
  };

  /// Construct the seeding algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  GridTripletSeedingAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

  /// Run the seeding algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  Acts::Experimental::CylindricalSpacePointGrid2::Config m_gridConfig;

  std::unique_ptr<const Acts::GridBinFinder<3ul>> m_bottomBinFinder{nullptr};
  std::unique_ptr<const Acts::GridBinFinder<3ul>> m_topBinFinder{nullptr};
  std::optional<Acts::Experimental::BroadTripletSeedFinder> m_seedFinder;
  std::optional<Acts::Experimental::BroadTripletSeedFilter> m_seedFilter;

  Acts::Delegate<bool(const SimSpacePoint&)> m_spacePointSelector{
      Acts::DelegateFuncTag<voidSpacePointSelector>{}};

  static bool voidSpacePointSelector(const SimSpacePoint& /*sp*/) {
    return true;
  }

  ReadDataHandle<SimSpacePointContainer> m_inputSpacePoints{this,
                                                            "InputSpacePoints"};

  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};

  /// Get the proper radius validity range given a middle space point candidate.
  /// In case the radius range changes according to the z-bin we need to
  /// retrieve the proper range. We can do this computation only once, since all
  /// the middle space point candidates belong to the same z-bin
  /// @param spM space point candidate to be used as middle SP in a seed
  /// @param rMiddleSPRange range object containing the minimum and maximum r for middle SP for a certain z bin
  std::pair<float, float> retrieveRadiusRangeForMiddle(
      const Acts::Experimental::ConstSpacePointProxy2& spM,
      const Acts::Range1D<float>& rMiddleSPRange) const;

  static inline bool itkFastTrackingCuts(float bottomRadius, float cotTheta) {
    static float rMin = 45;
    static float cotThetaMax = 1.5;

    if (bottomRadius < rMin &&
        (cotTheta > cotThetaMax || cotTheta < -cotThetaMax)) {
      return false;
    }
    return true;
  }

  static inline bool itkFastTrackingSPselect(const SimSpacePoint& sp) {
    // At small r we remove points beyond |z| > 200.
    float r = sp.r();
    float zabs = std::abs(sp.z());
    if (zabs > 200. && r < 45.) {
      return false;
    }

    // Remove space points beyond eta=4 if their z is larger than the max seed
    // z0 (150.)
    float cotTheta = 27.2899;  // corresponds to eta=4
    if ((zabs - 150.) > cotTheta * r) {
      return false;
    }
    return true;
  }
};

}  // namespace ActsExamples

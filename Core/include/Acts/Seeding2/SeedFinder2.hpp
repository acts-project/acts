// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData2/SeedContainer2.hpp"
#include "Acts/EventData2/SpacePointContainer2.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Seeding/SeedFinderUtils.hpp"
#include "Acts/Seeding2/SeedFilter2.hpp"
#include "Acts/Seeding2/detail/CandidatesForMiddleSp2.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace Acts {

class SeedFinder2 {
 private:
  enum class SpacePointCandidateType { eBottom, eTop };
  enum class MeasurementInfo { eDefault, eDetailed };

 public:
  struct DerivedConfig;
  struct DerivedOptions;

  struct Config {
    std::shared_ptr<SeedFilter2> seedFilter;

    // Seeding parameters used in the space-point grid creation and bin finding

    // Geometry settings + detector ROI
    // (r, z, phi) range for limiting location of all measurements and grid
    // creation
    float phiMin = -std::numbers::pi_v<float>;
    float phiMax = std::numbers::pi_v<float>;
    float zMin = -2800 * UnitConstants::mm;
    float zMax = 2800 * UnitConstants::mm;
    float rMax = 600 * UnitConstants::mm;
    /// WARNING: if rMin is smaller than impactMax, the bin size will be 2*pi,
    /// which will make seeding very slow!
    float rMin = 33 * UnitConstants::mm;

    /// Vector containing the z-bin edges for non equidistant binning in z
    std::vector<float> zBinEdges;

    /// Order of z bins to loop over when searching for SPs
    std::vector<std::size_t> zBinsCustomLooping;

    /// Radial bin size used in space-point grid
    float binSizeR = 1. * UnitConstants::mm;

    // Seeding parameters used to define the region of interest for middle
    // space-point

    /// Radial range for middle space-point
    /// The range can be defined manually with (rMinMiddle, rMaxMiddle). If
    /// useVariableMiddleSPRange is set to false and the vector rRangeMiddleSP
    /// is empty, we use (rMinMiddle, rMaxMiddle) to cut the middle space-points
    float rMinMiddle = 60.f * UnitConstants::mm;
    float rMaxMiddle = 120.f * UnitConstants::mm;
    /// If useVariableMiddleSPRange is set to false, the vector rRangeMiddleSP
    /// can be used to define a fixed r range for each z bin: {{rMin, rMax},
    /// ...}
    bool useVariableMiddleSPRange = false;
    /// Range defined in vector for each z bin
    std::vector<std::vector<float>> rRangeMiddleSP;
    /// If useVariableMiddleSPRange is true, the radial range will be calculated
    /// based on the maximum and minimum r values of the space-points in the
    /// event and a deltaR (deltaRMiddleMinSPRange, deltaRMiddleMaxSPRange)
    float deltaRMiddleMinSPRange = 10. * UnitConstants::mm;
    float deltaRMiddleMaxSPRange = 10. * UnitConstants::mm;

    // Seeding parameters used to define the cuts on space-point doublets

    /// Minimum radial distance between two doublet components (prefer
    /// deltaRMinTopSP and deltaRMinBottomSP to set separate values for outer
    /// and inner space-points)
    float deltaRMin = 5 * UnitConstants::mm;
    /// Maximum radial distance between two doublet components (prefer
    /// deltaRMaxTopSP and deltaRMacBottomSP to set separate values for outer
    /// and inner space-points)
    float deltaRMax = 270 * UnitConstants::mm;
    /// Minimum radial distance between middle-outer doublet components
    float deltaRMinTopSP = std::numeric_limits<float>::quiet_NaN();
    /// Maximum radial distance between middle-outer doublet components
    float deltaRMaxTopSP = std::numeric_limits<float>::quiet_NaN();
    /// Minimum radial distance between inner-middle doublet components
    float deltaRMinBottomSP = std::numeric_limits<float>::quiet_NaN();
    /// Maximum radial distance between inner-middle doublet components
    float deltaRMaxBottomSP = std::numeric_limits<float>::quiet_NaN();

    /// Maximum value of z-distance between space-points in doublet
    float deltaZMax =
        std::numeric_limits<float>::infinity() * UnitConstants::mm;

    /// Maximum allowed cotTheta between two space-points in doublet, used to
    /// check if forward angle is within bounds
    float cotThetaMax = 10.01788;  // equivalent to eta = 3 (pseudorapidity)

    /// Limiting location of collision region in z-axis used to check if doublet
    /// origin is within reasonable bounds
    float collisionRegionMin = -150 * UnitConstants::mm;
    float collisionRegionMax = +150 * UnitConstants::mm;

    /// Enable cut on the compatibility between interaction point and doublet,
    /// this is an useful approximation to speed up the seeding
    bool interactionPointCut = false;

    // Seeding parameters used to define the cuts on space-point triplets

    /// Minimum transverse momentum (pT) used to check the r-z slope
    /// compatibility of triplets with maximum multiple scattering effect
    /// (produced by the minimum allowed pT particle) + a certain uncertainty
    /// term. Check the documentation for more information
    /// https://acts.readthedocs.io/en/latest/core/reconstruction/pattern_recognition/seeding.html
    float minPt = 400. * UnitConstants::MeV;
    /// Number of sigmas of scattering angle to be considered in the minimum pT
    /// scattering term
    float sigmaScattering = 5;
    /// Term that accounts for the thickness of scattering medium in radiation
    /// lengths in the Lynch & Dahl correction to the Highland equation default
    /// is 5%
    /// TODO: necessary to make amount of material dependent on detector region?
    float radLengthPerSeed = 0.05;
    /// Maximum transverse momentum for scattering calculation
    float maxPtScattering = 10 * UnitConstants::GeV;
    /// Maximum value of impact parameter estimation of the seed candidates
    float impactMax = 20. * UnitConstants::mm;
    /// Parameter which can loosen the tolerance of the track seed to form a
    /// helix. This is useful for e.g. misaligned seeding.
    float helixCutTolerance = 1.;

    // Seeding parameters used for quality seed confirmation

    /// Enable quality seed confirmation, this is different than default seeding
    /// confirmation because it can also be defined for different (r, z) regions
    /// of the detector (e.g. forward or central region) by
    /// SeedConfirmationRange. Seeds are classified as "high-quality" seeds and
    /// normal quality seeds. Normal quality seeds are only selected if no other
    /// "high-quality" seeds has been found for that inner-middle doublet.
    bool seedConfirmation = false;
    /// Contains parameters for central seed confirmation
    SeedConfirmationRangeConfig centralSeedConfirmationRange;
    /// Contains parameters for forward seed confirmation
    SeedConfirmationRangeConfig forwardSeedConfirmationRange;
    /// Maximum number (minus one) of accepted seeds per middle space-point
    unsigned int maxSeedsPerSpM = 5;

    /// If seedConfirmation is true we classify seeds as "high-quality" seeds.
    /// Seeds that are not confirmed as "high-quality" are only selected if no
    /// other "high-quality" seed has been found for that inner-middle doublet
    /// Maximum number of normal seeds (not classified as "high-quality" seeds)
    /// in seed confirmation
    std::size_t maxSeedsPerSpMConf = std::numeric_limits<std::size_t>::max();
    /// Maximum number of "high-quality" seeds for each inner-middle SP-dublet
    /// in seed confirmation. If the limit is reached we check if there is a
    /// normal quality seed to be replaced
    std::size_t maxQualitySeedsPerSpMConf =
        std::numeric_limits<std::size_t>::max();

    // Other parameters

    /// Alignment uncertainties, used for uncertainties in the
    /// non-measurement-plane of the modules
    /// which otherwise would be 0
    /// will be added to spacepoint measurement uncertainties (and therefore
    /// also multiplied by sigmaError)
    /// FIXME: call align1 and align2
    float zAlign = 0 * UnitConstants::mm;
    float rAlign = 0 * UnitConstants::mm;
    /// used for measurement (+alignment) uncertainties.
    /// find seeds within 5sigma error ellipse
    float sigmaError = 5;

    /// only for Cuda plugin
    int maxBlockSize = 1024;
    int nTrplPerSpBLimit = 100;
    int nAvgTrplPerSpBLimit = 2;

    /// Delegates for accessors to detailed information on double measurement
    /// that produced the space point. This is mainly referring to space points
    /// produced when combining measurement from strips on back-to-back modules.
    /// Enables setting of the following delegates.
    bool useDetailedDoubleMeasurementInfo = false;

    /// Tolerance parameter used to check the compatibility of space-point
    /// coordinates in xyz. This is only used in a detector specific check for
    /// strip modules
    float toleranceParam = 1.1 * UnitConstants::mm;

    // Delegate to apply experiment specific cuts during doublet finding
    Delegate<bool(float /*bottomRadius*/, float /*cotTheta*/)> experimentCuts{
        DelegateFuncTag<&noopExperimentCuts>{}};

    DerivedConfig derive() const;

   private:
    static inline bool noopExperimentCuts(float /*bottomRadius*/,
                                          float /*cotTheta*/) {
      return true;
    }
  };

  struct DerivedConfig : public Config {
    float highland = 0;
    float maxScatteringAngle2 = 0;
  };

  struct Options {
    /// location of beam in x,y plane.
    /// used as offset for Space Points
    Vector2 beamPos{0 * UnitConstants::mm, 0 * UnitConstants::mm};
    /// field induction
    float bFieldInZ = 2.08 * UnitConstants::T;

    Range1D<float> rMiddleSpRange;

    DerivedOptions derive(const DerivedConfig& config) const;
  };

  struct DerivedOptions : public Options {
    float pTPerHelixRadius = std::numeric_limits<float>::quiet_NaN();
    float minHelixDiameter2 = std::numeric_limits<float>::quiet_NaN();
    float pT2perRadius = std::numeric_limits<float>::quiet_NaN();
    float sigmapT2perRadius = std::numeric_limits<float>::quiet_NaN();
    float multipleScattering2 = std::numeric_limits<float>::quiet_NaN();
  };

  struct State {
    std::vector<std::size_t> bottomSpGroupOffets;
    std::vector<std::size_t> topSpGroupOffets;

    std::vector<SpacePointIndex2> compatibleBottomSp;
    std::vector<SpacePointIndex2> compatibleTopSp;
    // contains parameters required to calculate circle with linear equation
    // ...for bottom-middle
    std::vector<LinCircle> linCirclesBottom;
    // ...for middle-top
    std::vector<LinCircle> linCirclesTop;

    std::vector<std::uint32_t> sortedBottoms;
    std::vector<std::uint32_t> sortedTops;

    std::vector<float> linCircleCotThetaBottom;
    std::vector<float> linCircleCotThetaTop;

    std::vector<SpacePointIndex2> topSpVec;
    std::vector<float> curvatures;
    std::vector<float> impactParameters;

    CandidatesForMiddleSp2 candidatesCollector;
    std::vector<TripletCandidate2> sortedCandidates;

    SeedFilter2::Options filterOptions;
    SeedFilter2::State filterState;
  };

  explicit SeedFinder2(const DerivedConfig& config,
                       std::unique_ptr<const Logger> logger = getDefaultLogger(
                           "SeedFinder2", Logging::Level::INFO));

  const DerivedConfig& config() const { return m_cfg; }

  void createSeeds(
      const DerivedOptions& options, State& state,
      const SpacePointContainer2& spacePoints,
      const SpacePointColumn2<float>& rColumn,
      const SpacePointColumn2<float>* varianceRColumn,
      const SpacePointColumn2<float>* varianceZColumn,
      const std::vector<std::vector<SpacePointIndex2>>& bottomSpGroups,
      const std::vector<SpacePointIndex2>& middleSpGroup,
      const std::vector<std::vector<SpacePointIndex2>>& topSpGroups,
      SeedContainer2& outputSeeds) const;

 private:
  struct DubletCuts {
    /// minimum allowed r-distance between dublet components
    float deltaRMin;
    /// maximum allowed r-distance between dublet components
    float deltaRMax;
    /// minus one over radius of middle SP
    float uIP;
    /// square of uIP
    float uIP2;
    /// ratio between middle SP x position and radius
    float cosPhiM;
    /// ratio between middle SP y position and radius
    float sinPhiM;
  };

  DubletCuts deriveDubletCuts(const ConstSpacePointProxy2& spM,
                              const SpacePointColumn2<float>& rColumn) const;

  std::pair<float, float> retrieveRadiusRangeForMiddle(
      const ConstSpacePointProxy2& spM,
      const Range1D<float>& rMiddleSPRange) const;

  template <SpacePointCandidateType candidate_type>
  void createCompatibleDoublets(
      const DerivedOptions& options, const DubletCuts& cuts,
      const SpacePointContainer2& spacePoints,
      const SpacePointColumn2<float>& rColumn,
      const SpacePointColumn2<float>* varianceRColumn,
      const SpacePointColumn2<float>* varianceZColumn,
      const ConstSpacePointProxy2& middleSp,
      const std::vector<std::vector<SpacePointIndex2>>& candidateSpGroups,
      std::vector<std::size_t>& candidateSpGroupOffsets,
      std::vector<SpacePointIndex2>& compatibleSp,
      std::vector<LinCircle>& linCircles) const;

  template <MeasurementInfo measurement_info>
  void filterCandidates(
      const DerivedOptions& options, State& state,
      const SpacePointContainer2& spacePoints,
      const SpacePointColumn2<float>& rColumn,
      const SpacePointColumn2<float>* varianceRColumn,
      const SpacePointColumn2<float>* varianceZColumn,
      const SpacePointColumn2<Acts::Vector3>* topStripVectorColumn,
      const SpacePointColumn2<Acts::Vector3>* bottomStripVectorColumn,
      const SpacePointColumn2<Acts::Vector3>* stripCenterDistanceColumn,
      const SpacePointColumn2<Acts::Vector3>* topStripCenterPositionColumn,
      const ConstSpacePointProxy2& spM) const;

  /// check the compatibility of SPs coordinates in xyz assuming the
  /// Bottom-Middle direction with the strip measurement details
  bool xyzCoordinateCheck(
      const ConstSpacePointProxy2& sp,
      const SpacePointColumn2<Acts::Vector3>& topStripVectorColumn,
      const SpacePointColumn2<Acts::Vector3>& bottomStripVectorColumn,
      const SpacePointColumn2<Acts::Vector3>& stripCenterDistanceColumn,
      const SpacePointColumn2<Acts::Vector3>& topStripCenterPositionColumn,
      const double* spacepointPosition, double* outputCoordinates) const;

 private:
  const Logger& logger() const { return *m_logger; }

  DerivedConfig m_cfg;
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts

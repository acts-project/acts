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
#include "Acts/Seeding2/TripletSeedFilter2.hpp"
#include "Acts/Seeding2/detail/CandidatesForMiddleSp2.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace Acts {

/// @brief Triplet seeding algorithm front-end
///
/// This class implements the triplet seeding algorithm, which is typical
/// procedure to find track seeds using space points in a cylindrical detector.
/// It is designed to be fast, flexible, and configurable.
///
/// The algorithm works by first finding compatible doublets of space points,
/// two space points that can be connected by a track coming from the
/// interaction region, and then forming triplets by combinding these
/// doublets at a common middle space point. The triplets are then filtered
/// using a seed filter to produce a set of track seeds.
///
/// Note that this algorithm is designed and tuned for cylindrical detectors and
/// uses R-Z coordinates for the space points.
class TripletSeedFinder2 {
 private:
  enum class SpacePointCandidateType { eBottom, eTop };
  enum class MeasurementInfo { eDefault, eDetailed };

 public:
  struct DerivedConfig;
  struct DerivedOptions;

  struct Config {
    std::shared_ptr<TripletSeedFilter2> seedFilter;

    // Seeding parameters used in the space-point grid creation and bin finding

    /// Vector containing the z-bin edges for non equidistant binning in z
    std::vector<float> zBinEdges;

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

    // Seeding parameters used to define the cuts on space-point doublets

    /// Minimum radial distance between middle-outer doublet components
    float deltaRMinTopSP = 5 * UnitConstants::mm;
    /// Maximum radial distance between middle-outer doublet components
    float deltaRMaxTopSP = 270 * UnitConstants::mm;
    /// Minimum radial distance between inner-middle doublet components
    float deltaRMinBottomSP = 5 * UnitConstants::mm;
    /// Maximum radial distance between inner-middle doublet components
    float deltaRMaxBottomSP = 270 * UnitConstants::mm;

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

    /// Delegates for accessors to detailed information on double measurement
    /// that produced the space point. This is mainly referring to space points
    /// produced when combining measurement from strips on back-to-back modules.
    /// Enables setting of the following delegates.
    bool useDetailedDoubleMeasurementInfo = false;

    /// Tolerance parameter used to check the compatibility of space-point
    /// coordinates in xyz. This is only used in a detector specific check for
    /// strip modules
    float toleranceParam = 1.1 * UnitConstants::mm;

    /// Delegate to apply experiment specific cuts during doublet finding
    Delegate<bool(float /*bottomRadius*/, float /*cotTheta*/)> experimentCuts{
        DelegateFuncTag<&noopExperimentCuts>{}};

    /// Create a derived configuration object from the current configuration.
    /// Computes derived quantities and changes to internal seeding units.
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

    TripletSeedFilter2::State filterState;
  };

  explicit TripletSeedFinder2(const DerivedConfig& config,
                              std::unique_ptr<const Logger> logger =
                                  getDefaultLogger("TripletSeedFinder2",
                                                   Logging::Level::INFO));

  const DerivedConfig& config() const { return m_cfg; }

  /// Create all possible seeds from bottom, middle, and top space points.
  /// This function is the main entry point for the seeding algorithm.
  ///
  /// @param options frequently changing configuration (like beam position)
  /// @param state State object to cache memory and results
  /// @param spacePoints The container with space points
  /// @param rColumn The column with r values of space points
  /// @param varianceRColumn Optional column with r variance values
  /// @param varianceZColumn Optional column with z variance values
  /// @param bottomSps Group of space points to be used as innermost SP in a seed.
  /// @param middleSps Group of space points to be used as middle SP in a seed.
  /// @param topSps Group of space points to be used as outermost SP in a seed.
  /// @param outputSeeds Output container for the seeds
  void createSeeds(
      const DerivedOptions& options, State& state,
      const SpacePointContainer2& spacePoints,
      const SpacePointContainer2::DenseColumn<float>& rColumn,
      const SpacePointContainer2::DenseColumn<float>* varianceRColumn,
      const SpacePointContainer2::DenseColumn<float>* varianceZColumn,
      std::vector<SpacePointIndex2>& bottomSps,
      std::vector<SpacePointIndex2>& middleSps,
      std::vector<SpacePointIndex2>& topSps, SeedContainer2& outputSeeds) const;

 private:
  struct DoubletCuts {
    /// minimum allowed r-distance between doublet components
    float deltaRMin{};
    /// maximum allowed r-distance between doublet components
    float deltaRMax{};
    /// minus one over radius of middle SP
    float uIP{};
    /// square of uIP
    float uIP2{};
    /// ratio between middle SP x position and radius
    float cosPhiM{};
    /// ratio between middle SP y position and radius
    float sinPhiM{};
  };

  void deriveDoubletCuts(
      DoubletCuts& cuts, const ConstSpacePointProxy2& spM,
      const SpacePointContainer2::DenseColumn<float>& rColumn) const;

  /// Get the proper radius validity range given a middle space point candidate.
  /// In case the radius range changes according to the z-bin we need to
  /// retrieve the proper range. We can do this computation only once, since all
  /// the middle space point candidates belong to the same z-bin
  /// @param spM space point candidate to be used as middle SP in a seed
  /// @param rMiddleSPRange range object containing the minimum and maximum r for middle SP for a certain z bin.
  std::pair<float, float> retrieveRadiusRangeForMiddle(
      const ConstSpacePointProxy2& spM,
      const Range1D<float>& rMiddleSPRange) const;

  /// Iterates over dublets and tests the compatibility by applying a series of
  /// cuts that can be tested with only two SPs.
  ///
  /// @tparam candidate_type Type of space point candidate (e.g. Bottom or Top)
  ///
  /// @param config the configuration for the SeedFinder
  /// @param options frequently changing configuration (like beam position)
  /// @param cuts Doublet cuts that define the compatibility of space points
  /// @param spacePoints The container with space points
  /// @param rColumn The column with r values of space points
  /// @param varianceRColumn Optional column with r variance values
  /// @param varianceZColumn Optional column with z variance values
  /// @param middleSp Space point candidate to be used as middle SP in a seed
  /// @param candidateSps Group of space points to be used as candidates for
  /// middle SP in a seed.
  /// @param compatibleSp Output vector of compatible space points
  /// @param linCircles Output vector of LinCircle objects for the bottom-middle
  /// doublets
  /// @param cotThetas Output vector of cotTheta values for the bottom-middle
  /// doublets
  template <SpacePointCandidateType candidate_type>
  static void createDoublets(
      const DerivedConfig& config, const DerivedOptions& options,
      const DoubletCuts& cuts, const SpacePointContainer2& spacePoints,
      const SpacePointContainer2::DenseColumn<float>& rColumn,
      const SpacePointContainer2::DenseColumn<float>* varianceRColumn,
      const SpacePointContainer2::DenseColumn<float>* varianceZColumn,
      const ConstSpacePointProxy2& middleSp,
      const std::vector<SpacePointIndex2>& candidateSps,
      std::vector<SpacePointIndex2>& compatibleSp,
      std::vector<LinCircle>& linCircles, std::vector<float>& cotThetas);

  /// Create triplets from the bottom, middle, and top space points.
  ///
  /// @tparam measurement_info Type of measurement information (e.g. Default or Detailed)
  ///
  /// @param config the configuration for the SeedFinder
  /// @param options frequently changing configuration (like beam position)
  /// @param filterOptions Options for the seed filter
  /// @param filterState State object that holds memory used in SeedFilter
  /// @param spacePoints The container with space points
  /// @param rColumn The column with r values of space points
  /// @param varianceRColumn Optional column with r variance values
  /// @param varianceZColumn Optional column with z variance values
  /// @param topStripVectorColumn Optional column with top strip vectors
  /// @param bottomStripVectorColumn Optional column with bottom strip vectors
  /// @param stripCenterDistanceColumn Optional column with strip center distances
  /// @param topStripCenterPositionColumn Optional column with top strip center positions
  /// @param spM Space point candidate to be used as middle SP in a seed
  /// @param bottomSps Group of space points to be used as innermost SP in a seed.
  /// @param bottomLinCircles Vector containing bottom-middle SP parameters after
  /// reference frame transformation to the u-v space
  /// @param bottomCotThetas Vector containing cotTheta values for the bottom-middle SPs
  /// @param topSps Group of space points to be used as outermost SP in a seed.
  /// @param topLinCircles Vector containing middle-top SP parameters after
  /// reference frame transformation to the u-v space
  /// @param topCotThetas Vector containing cotTheta values for the middle-top SPs
  /// @param sortedBottoms Output vector of sorted bottom space points
  /// @param sortedTops Output vector of sorted top space points
  /// @param topSpVec Output vector of top space point indices
  /// @param curvatures Output vector of curvature values for the triplets
  /// @param impactParameters Output vector of impact parameter values for the triplets
  /// @param candidatesCollector Collector for candidates for middle space points
  template <MeasurementInfo measurement_info>
  static void createTriplets(
      const DerivedConfig& config, const DerivedOptions& options,
      const TripletSeedFilter2::Options& filterOptions,
      TripletSeedFilter2::State& filterState,
      const SpacePointContainer2& spacePoints,
      const SpacePointContainer2::DenseColumn<float>& rColumn,
      const SpacePointContainer2::DenseColumn<float>* varianceRColumn,
      const SpacePointContainer2::DenseColumn<float>* varianceZColumn,
      const SpacePointContainer2::DenseColumn<Vector3>* topStripVectorColumn,
      const SpacePointContainer2::DenseColumn<Vector3>* bottomStripVectorColumn,
      const SpacePointContainer2::DenseColumn<Vector3>*
          stripCenterDistanceColumn,
      const SpacePointContainer2::DenseColumn<Vector3>*
          topStripCenterPositionColumn,
      const ConstSpacePointProxy2& spM,
      const std::vector<SpacePointIndex2>& bottomSps,
      const std::vector<LinCircle>& bottomLinCircles,
      const std::vector<float>& bottomCotThetas,
      const std::vector<SpacePointIndex2>& topSps,
      const std::vector<LinCircle>& topLinCircles,
      const std::vector<float>& topCotThetas,
      std::vector<std::uint32_t>& sortedBottoms,
      std::vector<std::uint32_t>& sortedTops,
      std::vector<SpacePointIndex2>& topSpVec, std::vector<float>& curvatures,
      std::vector<float>& impactParameters,
      CandidatesForMiddleSp2& candidatesCollector);

  /// Check the compatibility of SPs coordinates in xyz assuming the
  /// Bottom-Middle direction with the strip measurement details
  static bool xyzCoordinateCheck(
      double toleranceParam, const ConstSpacePointProxy2& sp,
      const SpacePointContainer2::DenseColumn<Vector3>& topStripVectorColumn,
      const SpacePointContainer2::DenseColumn<Vector3>& bottomStripVectorColumn,
      const SpacePointContainer2::DenseColumn<Vector3>&
          stripCenterDistanceColumn,
      const SpacePointContainer2::DenseColumn<Vector3>&
          topStripCenterPositionColumn,
      const double* spacepointPosition, double* outputCoordinates);

 private:
  const Logger& logger() const { return *m_logger; }

  DerivedConfig m_cfg;
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts

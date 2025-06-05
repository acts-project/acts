// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
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
 public:
  enum class SpacePointCandidateType { eBottom, eTop };

  struct DerivedConfig;
  struct DerivedOptions;

  struct Config {
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

    /// Delegate to apply experiment specific cuts during doublet finding
    Delegate<bool(float /*bottomRadius*/, float /*cotTheta*/)> experimentCuts;

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

    std::shared_ptr<TripletSeedFilter2> filter;

    /// Create a derived configuration object from the current configuration.
    /// Computes derived quantities and changes to internal seeding units.
    DerivedConfig derive() const;
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

  struct DoubletCuts {
    /// Maximum value of impact parameter estimation of the seed candidates
    float impactMax = 20. * UnitConstants::mm;

    /// Enable cut on the compatibility between interaction point and doublet,
    /// this is an useful approximation to speed up the seeding
    bool interactionPointCut = false;

    /// minimum allowed r-distance between doublet components
    float deltaRMin = 5 * UnitConstants::mm;
    /// maximum allowed r-distance between doublet components
    float deltaRMax = 270 * UnitConstants::mm;

    /// Limiting location of collision region in z-axis used to check if doublet
    /// origin is within reasonable bounds
    float collisionRegionMin = -150 * UnitConstants::mm;
    float collisionRegionMax = +150 * UnitConstants::mm;

    /// Maximum allowed cotTheta between two space-points in doublet, used to
    /// check if forward angle is within bounds
    float cotThetaMax = 10.01788;  // equivalent to eta = 3 (pseudorapidity)

    /// Maximum value of z-distance between space-points in doublet
    float deltaZMax =
        std::numeric_limits<float>::infinity() * UnitConstants::mm;

    float minHelixDiameter2 = std::numeric_limits<float>::quiet_NaN();

    /// Delegate to apply experiment specific cuts during doublet finding
    Delegate<bool(float /*bottomRadius*/, float /*cotTheta*/)> experimentCuts;
  };

  struct TripletCuts {
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

    float highland = 0;
    float pTPerHelixRadius = std::numeric_limits<float>::quiet_NaN();
    float minHelixDiameter2 = std::numeric_limits<float>::quiet_NaN();
    float sigmapT2perRadius = std::numeric_limits<float>::quiet_NaN();
    float multipleScattering2 = std::numeric_limits<float>::quiet_NaN();

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

    /// Tolerance parameter used to check the compatibility of space-point
    /// coordinates in xyz. This is only used in a detector specific check for
    /// strip modules
    float toleranceParam = 1.1 * UnitConstants::mm;
  };

  struct State {
    Options options;

    DoubletCuts bottomDoubletCuts;
    DoubletCuts topDoubletCuts;

    TripletCuts tripletCuts;

    TripletSeedFilter2::State filter;
  };

  struct Doublets {
    std::vector<SpacePointIndex2> spacePoints;
    /// contains parameters required to calculate a circle with linear equation
    std::vector<LinCircle> linCircles;
    std::vector<float> cotTheta;

    [[nodiscard]] bool empty() const { return spacePoints.empty(); }

    [[nodiscard]] std::size_t size() const { return spacePoints.size(); }

    void clear() {
      spacePoints.clear();
      linCircles.clear();
      cotTheta.clear();
    }

    void emplace_back(SpacePointIndex2 sp, const LinCircle& linCircle) {
      spacePoints.emplace_back(sp);
      linCircles.emplace_back(linCircle);
      cotTheta.emplace_back(linCircle.cotTheta);
    }
  };

  struct MiddleSpacePointInfo {
    /// minus one over radius of middle SP
    float uIP{};
    /// square of uIP
    float uIP2{};
    /// ratio between middle SP x position and radius
    float cosPhiM{};
    /// ratio between middle SP y position and radius
    float sinPhiM{};
  };

  struct TripletCache {
    std::vector<std::uint32_t> sortedBottoms;
    std::vector<std::uint32_t> sortedTops;
  };

  struct TripletTopCandidates {
    std::vector<SpacePointIndex2> topSpacePoints;
    std::vector<float> curvatures;
    std::vector<float> impactParameters;

    void resize(std::size_t size) {
      topSpacePoints.resize(size);
      curvatures.resize(size);
      impactParameters.resize(size);
    }

    void clear() {
      topSpacePoints.clear();
      curvatures.clear();
      impactParameters.clear();
    }

    void emplace_back(SpacePointIndex2 spT, float curvature,
                      float impactParameter) {
      topSpacePoints.emplace_back(spT);
      curvatures.emplace_back(curvature);
      impactParameters.emplace_back(impactParameter);
    }
  };

  struct Cache {
    TripletSeedFilter2::Cache filter;

    Doublets bottomDoublets;
    Doublets topDoublets;

    TripletCache tripletCache;
    TripletTopCandidates tripletTopCandidates;

    CandidatesForMiddleSp2 candidatesCollector;
    std::vector<TripletCandidate2> sortedCandidates;
  };

  /// Collection of pointers to the space point container and its
  /// additional columns. This is used as a basket to pass around
  /// the input data for the triplet seed finder.
  class ContainerPointers {
   public:
    /// Minimal input: space points and r column.
    ContainerPointers(const SpacePointContainer2& spacePoints,
                      const SpacePointContainer2::DenseColumn<float>& rColumn)
        : m_spacePoints(&spacePoints), m_rColumn(&rColumn) {}

    /// Space points, r column, and variance columns.
    ContainerPointers(
        const SpacePointContainer2& spacePoints,
        const SpacePointContainer2::DenseColumn<float>& rColumn,
        const SpacePointContainer2::DenseColumn<float>& varianceRColumn,
        const SpacePointContainer2::DenseColumn<float>& varianceZColumn)
        : m_spacePoints(&spacePoints),
          m_rColumn(&rColumn),
          m_varianceRColumn(&varianceRColumn),
          m_varianceZColumn(&varianceZColumn) {}

    /// Space points, r column, variance columns, and strip columns.
    ContainerPointers(
        const SpacePointContainer2& spacePoints,
        const SpacePointContainer2::DenseColumn<float>& rColumn,
        const SpacePointContainer2::DenseColumn<float>& varianceRColumn,
        const SpacePointContainer2::DenseColumn<float>& varianceZColumn,
        const SpacePointContainer2::DenseColumn<Vector3>& topStripVectorColumn,
        const SpacePointContainer2::DenseColumn<Vector3>&
            bottomStripVectorColumn,
        const SpacePointContainer2::DenseColumn<Vector3>&
            stripCenterDistanceColumn,
        const SpacePointContainer2::DenseColumn<Vector3>&
            topStripCenterPositionColumn)
        : m_spacePoints(&spacePoints),
          m_rColumn(&rColumn),
          m_varianceRColumn(&varianceRColumn),
          m_varianceZColumn(&varianceZColumn),
          m_topStripVectorColumn(&topStripVectorColumn),
          m_bottomStripVectorColumn(&bottomStripVectorColumn),
          m_stripCenterDistanceColumn(&stripCenterDistanceColumn),
          m_topStripCenterPositionColumn(&topStripCenterPositionColumn) {}

    /// Pointer to the copied-from index column, if available.
    const SpacePointContainer2::DenseColumn<SpacePointIndex2>*
        copiedFromIndexColumn = nullptr;

    [[nodiscard]] const SpacePointContainer2& spacePoints() const {
      return *m_spacePoints;
    }
    [[nodiscard]] const SpacePointContainer2::DenseColumn<float>& rColumn()
        const {
      return *m_rColumn;
    }
    [[nodiscard]] const SpacePointContainer2::DenseColumn<float>&
    varianceRColumn() const {
      return *m_varianceRColumn;
    }
    [[nodiscard]] const SpacePointContainer2::DenseColumn<float>&
    varianceZColumn() const {
      return *m_varianceZColumn;
    }
    [[nodiscard]] const SpacePointContainer2::DenseColumn<Vector3>&
    topStripVectorColumn() const {
      return *m_topStripVectorColumn;
    }
    [[nodiscard]] const SpacePointContainer2::DenseColumn<Vector3>&
    bottomStripVectorColumn() const {
      return *m_bottomStripVectorColumn;
    }
    [[nodiscard]] const SpacePointContainer2::DenseColumn<Vector3>&
    stripCenterDistanceColumn() const {
      return *m_stripCenterDistanceColumn;
    }
    [[nodiscard]] const SpacePointContainer2::DenseColumn<Vector3>&
    topStripCenterPositionColumn() const {
      return *m_topStripCenterPositionColumn;
    }

    [[nodiscard]] bool hasVarianceColumns() const {
      return m_varianceRColumn != nullptr && m_varianceZColumn != nullptr;
    }
    [[nodiscard]] bool hasStripColumns() const {
      return m_topStripVectorColumn != nullptr &&
             m_bottomStripVectorColumn != nullptr &&
             m_stripCenterDistanceColumn != nullptr &&
             m_topStripCenterPositionColumn != nullptr;
    }
    [[nodiscard]] bool hasCopiedFromIndexColumn() const {
      return copiedFromIndexColumn != nullptr;
    }

   private:
    const SpacePointContainer2* m_spacePoints = nullptr;
    const SpacePointContainer2::DenseColumn<float>* m_rColumn = nullptr;

    const SpacePointContainer2::DenseColumn<float>* m_varianceRColumn = nullptr;
    const SpacePointContainer2::DenseColumn<float>* m_varianceZColumn = nullptr;

    const SpacePointContainer2::DenseColumn<Vector3>* m_topStripVectorColumn =
        nullptr;
    const SpacePointContainer2::DenseColumn<Vector3>*
        m_bottomStripVectorColumn = nullptr;
    const SpacePointContainer2::DenseColumn<Vector3>*
        m_stripCenterDistanceColumn = nullptr;
    const SpacePointContainer2::DenseColumn<Vector3>*
        m_topStripCenterPositionColumn = nullptr;
  };

  explicit TripletSeedFinder2(const DerivedConfig& config,
                              std::unique_ptr<const Logger> logger =
                                  getDefaultLogger("TripletSeedFinder2",
                                                   Logging::Level::INFO));

  const DerivedConfig& config() const { return m_cfg; }

  /// Initialize the state of the seed finder with the provided options.
  ///
  /// @param state State of the seed finder
  /// @param options Derived options that may include configuration parameters
  /// @note This function should be called before using the seed finder
  void initialize(State& state, const DerivedOptions& options) const;

  /// Create all possible seeds from bottom, middle, and top space points.
  /// This function is the main entry point for the seeding algorithm.
  ///
  /// @param options frequently changing configuration (like beam position)
  /// @param state State of the seed finder
  /// @param cache Cache object to store intermediate results
  /// @param containerPointers Space point container and its extra columns
  /// @param bottomSps Group of space points to be used as innermost SP in a seed
  /// @param middleSps Group of space points to be used as middle SP in a seed
  /// @param topSps Group of space points to be used as outermost SP in a seed
  /// @param outputSeeds Output container for the seeds
  void createSeeds(State& state, Cache& cache,
                   const ContainerPointers& containerPointers,
                   std::span<const SpacePointIndex2> bottomSps,
                   std::span<const SpacePointIndex2> middleSps,
                   std::span<const SpacePointIndex2> topSps,
                   SeedContainer2& outputSeeds) const;

  /// Create all possible seeds from bottom, middle, and top space points.
  ///
  /// @param options frequently changing configuration (like beam position)
  /// @param state State of the seed finder
  /// @param cache Cache object to store intermediate results
  /// @param containerPointers Space point container and its extra columns
  /// @param bottomSps Group of space points to be used as innermost SP in a seed
  /// @param middleSp Space point candidate to be used as middle SP in a seed
  /// @param topSps Group of space points to be used as outermost SP in a seed
  /// @param outputSeeds Output container for the seeds
  void createSeeds(State& state, Cache& cache,
                   const ContainerPointers& containerPointers,
                   std::span<const SpacePointIndex2> bottomSps,
                   SpacePointIndex2 middleSp,
                   std::span<const SpacePointIndex2> topSps,
                   SeedContainer2& outputSeeds) const;

  /// Derives doublet cuts based on the provided options and space point
  /// candidate type (bottom or top).
  ///
  /// @param options Derived options that may include configuration parameters
  /// @param spacePointCandidateType Type of space point candidate (e.g. Bottom or Top)
  /// @return A DoubletCuts object containing the derived cuts
  DoubletCuts deriveDoubletCuts(
      const DerivedOptions& options,
      SpacePointCandidateType spacePointCandidateType) const;

  /// Derives triplet cuts based on the provided options.
  ///
  /// @param options Derived options that may include configuration parameters
  /// @return A TripletCuts object containing the derived cuts
  TripletCuts deriveTripletCuts(const DerivedOptions& options) const;

  /// Iterates over dublets and tests the compatibility by applying a series of
  /// cuts that can be tested with only two SPs.
  ///
  /// @tparam candidate_type Type of space point candidate (e.g. Bottom or Top)
  ///
  /// @param cuts Doublet cuts that define the compatibility of space points
  /// @param containerPointers Space point container and its extra columns
  /// @param middleSp Space point candidate to be used as middle SP in a seed
  /// @param middleSpInfo Information about the middle space point
  /// @param candidateSps Group of space points to be used as candidates for
  ///                     middle SP in a seed
  /// @param compatibleDoublets Output container for compatible doublets
  template <SpacePointCandidateType candidate_type>
  static void createDoublets(const DoubletCuts& cuts,
                             const ContainerPointers& containerPointers,
                             const ConstSpacePointProxy2& middleSp,
                             const MiddleSpacePointInfo& middleSpInfo,
                             std::span<const SpacePointIndex2> candidateSps,
                             Doublets& compatibleDoublets);

  /// Create triplets from the bottom, middle, and top space points.
  ///
  /// @tparam measurement_info Type of measurement information (e.g. Default or Detailed)
  ///
  /// @param cache Cache object to store intermediate results
  /// @param cuts Triplet cuts that define the compatibility of space points
  /// @param filter Triplet seed filter that defines the filtering criteria
  /// @param filterState State object that holds the state of the filter
  /// @param filterCache Cache object that holds memory used in SeedFilter
  /// @param containerPointers Space point container and its extra columns
  /// @param spM Space point candidate to be used as middle SP in a seed
  /// @param bottomDoublets Bottom doublets to be used for triplet creation
  /// @param topDoublets Top doublets to be used for triplet creation
  /// @param tripletTopCandidates Cache for triplet top candidates
  /// @param candidatesCollector Collector for candidates for middle space points
  static void createTriplets(TripletCache& cache, const TripletCuts& cuts,
                             const TripletSeedFilter2& filter,
                             const TripletSeedFilter2::Options& filterOptions,
                             TripletSeedFilter2::State& filterState,
                             TripletSeedFilter2::Cache& filterCache,
                             const ContainerPointers& containerPointers,
                             const ConstSpacePointProxy2& spM,
                             const Doublets& bottomDoublets,
                             const Doublets& topDoublets,
                             TripletTopCandidates& tripletTopCandidates,
                             CandidatesForMiddleSp2& candidatesCollector);

  /// Create triplets from the bottom, middle, and top space points.
  ///
  /// @tparam measurement_info Type of measurement information (e.g. Default or Detailed)
  ///
  /// @param cache Cache object to store intermediate results
  /// @param cuts Triplet cuts that define the compatibility of space points
  /// @param filter Triplet seed filter that defines the filtering criteria
  /// @param filterState State object that holds the state of the filter
  /// @param filterCache Cache object that holds memory used in SeedFilter
  /// @param containerPointers Space point container and its extra columns
  /// @param spM Space point candidate to be used as middle SP in a seed
  /// @param bottomDoublets Bottom doublets to be used for triplet creation
  /// @param topDoublets Top doublets to be used for triplet creation
  /// @param tripletTopCandidates Cache for triplet top candidates
  /// @param candidatesCollector Collector for candidates for middle space points
  static void createTripletsDetailed(
      const TripletCuts& cuts, const TripletSeedFilter2& filter,
      const TripletSeedFilter2::Options& filterOptions,
      TripletSeedFilter2::State& filterState,
      TripletSeedFilter2::Cache& filterCache,
      const ContainerPointers& containerPointers,
      const ConstSpacePointProxy2& spM, const Doublets& bottomDoublets,
      const Doublets& topDoublets, TripletTopCandidates& tripletTopCandidates,
      CandidatesForMiddleSp2& candidatesCollector);

 private:
  static MiddleSpacePointInfo computeMiddleSpacePointInfo(
      const ConstSpacePointProxy2& spM,
      const SpacePointContainer2::DenseColumn<float>& rColumn);

  /// Get the proper radius validity range given a middle space point candidate.
  /// In case the radius range changes according to the z-bin we need to
  /// retrieve the proper range. We can do this computation only once, since all
  /// the middle space point candidates belong to the same z-bin
  /// @param spM space point candidate to be used as middle SP in a seed
  /// @param rMiddleSPRange range object containing the minimum and maximum r for middle SP for a certain z bin
  std::pair<float, float> retrieveRadiusRangeForMiddle(
      const ConstSpacePointProxy2& spM,
      const Range1D<float>& rMiddleSPRange) const;

  /// Check the compatibility of strip space point coordinates in xyz assuming
  /// the Bottom-Middle direction with the strip measurement details
  static bool stripCoordinateCheck(float tolerance,
                                   const ConstSpacePointProxy2& sp,
                                   const ContainerPointers& containerPointers,
                                   const Vector3& spacePointPosition,
                                   Vector3& outputCoordinates);

 private:
  const Logger& logger() const { return *m_logger; }

  DerivedConfig m_cfg;
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Acts

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SpacePointMutableData.hpp"
#include "Acts/Seeding/CandidatesForMiddleSp.hpp"
#include "Acts/Seeding/Neighbour.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinderUtils.hpp"
#include "Acts/Seeding/detail/UtilityFunctions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"

#include <memory>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

namespace Acts {

template <typename Coll>
concept GridBinCollection =
    std::ranges::random_access_range<Coll> &&
    std::same_as<typename Coll::value_type, std::size_t>;

template <typename collection_t, typename external_t, std::size_t N = 3ul>
concept CollectionStoresSeedsTo =
    requires(collection_t coll, Seed<external_t, N> seed) {
      detail::pushBackOrInsertAtEnd(coll, seed);
    };

enum class SpacePointCandidateType : short { eBottom, eTop };

enum class DetectorMeasurementInfo : short { eDefault, eDetailed };

/// @brief Structure that holds configuration parameters for the seed finder algorithm
template <typename SpacePoint>
struct SeedFinderConfig {
  std::shared_ptr<SeedFilter<SpacePoint>> seedFilter;

  /// Seeding parameters used in the space-point grid creation and bin finding

  /// Geometry Settings + Detector ROI
  /// (r, z, phi) range for limiting location of all measurements and grid
  /// creation
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
  std::vector<std::size_t> zBinsCustomLooping = {};

  /// Radial bin size used in space-point grid
  float binSizeR = 1. * UnitConstants::mm;

  /// Seeding parameters used to define the region of interest for middle
  /// space-point

  /// Radial range for middle space-point
  /// The range can be defined manually with (rMinMiddle, rMaxMiddle). If
  /// useVariableMiddleSPRange is set to false and the vector rRangeMiddleSP is
  /// empty, we use (rMinMiddle, rMaxMiddle) to cut the middle space-points
  float rMinMiddle = 60.f * UnitConstants::mm;
  float rMaxMiddle = 120.f * UnitConstants::mm;
  /// If useVariableMiddleSPRange is set to false, the vector rRangeMiddleSP can
  /// be used to define a fixed r range for each z bin: {{rMin, rMax}, ...}
  bool useVariableMiddleSPRange = false;
  /// Range defined in vector for each z bin
  std::vector<std::vector<float>> rRangeMiddleSP;
  /// If useVariableMiddleSPRange is true, the radial range will be calculated
  /// based on the maximum and minimum r values of the space-points in the event
  /// and a deltaR (deltaRMiddleMinSPRange, deltaRMiddleMaxSPRange)
  float deltaRMiddleMinSPRange = 10. * UnitConstants::mm;
  float deltaRMiddleMaxSPRange = 10. * UnitConstants::mm;

  /// Seeding parameters used to define the cuts on space-point doublets

  /// Minimum radial distance between two doublet components (prefer
  /// deltaRMinTopSP and deltaRMinBottomSP to set separate values for outer and
  /// inner space-points)
  float deltaRMin = 5 * UnitConstants::mm;
  /// Maximum radial distance between two doublet components (prefer
  /// deltaRMaxTopSP and deltaRMacBottomSP to set separate values for outer and
  /// inner space-points)
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
  float deltaZMax = std::numeric_limits<float>::infinity() * UnitConstants::mm;

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

  /// Seeding parameters used to define the cuts on space-point triplets

  /// Minimum transverse momentum (pT) used to check the r-z slope compatibility
  /// of triplets with maximum multiple scattering effect (produced by the
  /// minimum allowed pT particle) + a certain uncertainty term. Check the
  /// documentation for more information
  /// https://acts.readthedocs.io/en/latest/core/reconstruction/pattern_recognition/seeding.html
  float minPt = 400. * UnitConstants::MeV;
  /// Number of sigmas of scattering angle to be considered in the minimum pT
  /// scattering term
  float sigmaScattering = 5;
  /// Term that accounts for the thickness of scattering medium in radiation
  /// lengths in the Lynch & Dahl correction to the Highland equation default is
  /// 5%
  /// TODO: necessary to make amount of material dependent on detector region?
  float radLengthPerSeed = 0.05;
  /// Maximum transverse momentum for scattering calculation
  float maxPtScattering = 10 * UnitConstants::GeV;
  /// Maximum value of impact parameter estimation of the seed candidates
  float impactMax = 20. * UnitConstants::mm;
  /// Parameter which can loosen the tolerance of the track seed to form a
  /// helix. This is useful for e.g. misaligned seeding.
  float helixCutTolerance = 1.;

  /// Seeding parameters used for quality seed confirmation

  /// Enable quality seed confirmation, this is different than default seeding
  /// confirmation because it can also be defined for different (r, z) regions
  /// of the detector (e.g. forward or central region) by SeedConfirmationRange.
  /// Seeds are classified as "high-quality" seeds and normal quality seeds.
  /// Normal quality seeds are only selected if no other "high-quality" seeds
  /// has been found for that inner-middle doublet.
  bool seedConfirmation = false;
  /// Contains parameters for central seed confirmation
  SeedConfirmationRangeConfig centralSeedConfirmationRange;
  /// Contains parameters for forward seed confirmation
  SeedConfirmationRangeConfig forwardSeedConfirmationRange;
  /// Maximum number (minus one) of accepted seeds per middle space-point
  unsigned int maxSeedsPerSpM = 5;

  /// Other parameters

  /// Alignment uncertainties, used for uncertainties in the
  /// non-measurement-plane of the modules
  /// which otherwise would be 0
  /// will be added to spacepoint measurement uncertainties (and therefore also
  /// multiplied by sigmaError)
  /// FIXME: call align1 and align2
  float zAlign = 0 * UnitConstants::mm;
  float rAlign = 0 * UnitConstants::mm;
  /// used for measurement (+alignment) uncertainties.
  /// find seeds within 5sigma error ellipse
  float sigmaError = 5;

  /// derived values, set on SeedFinder construction
  float highland = 0;
  float maxScatteringAngle2 = 0;

  /// only for Cuda plugin
  int maxBlockSize = 1024;
  int nTrplPerSpBLimit = 100;
  int nAvgTrplPerSpBLimit = 2;

  /// Delegates for accessors to detailed information on double measurement that
  /// produced the space point.
  /// This is mainly referring to space points produced when combining
  /// measurement from strips on back-to-back modules.
  /// Enables setting of the following delegates.
  bool useDetailedDoubleMeasurementInfo = false;

  Delegate<bool(const SpacePoint&)> spacePointSelector{
      DelegateFuncTag<voidSpacePointSelector>{}};

  static bool voidSpacePointSelector(const SpacePoint& /*sp*/) { return true; }

  /// Tolerance parameter used to check the compatibility of space-point
  /// coordinates in xyz. This is only used in a detector specific check for
  /// strip modules
  float toleranceParam = 1.1 * UnitConstants::mm;

  // Delegate to apply experiment specific cuts during doublet finding
  Delegate<bool(float /*bottomRadius*/, float /*cotTheta*/)> experimentCuts{
      DelegateFuncTag<&noopExperimentCuts>{}};

  bool isInInternalUnits = false;

  SeedFinderConfig toInternalUnits() const {
    if (isInInternalUnits) {
      throw std::runtime_error(
          "Repeated conversion to internal units for SeedFinderConfig");
    }
    // Make sure the shared ptr to the seed filter is not a nullptr
    // And make sure the seed filter config is in internal units as well
    if (!seedFilter) {
      throw std::runtime_error(
          "Invalid values for the seed filter inside the seed filter config: "
          "nullptr");
    }
    if (!seedFilter->getSeedFilterConfig().isInInternalUnits) {
      throw std::runtime_error(
          "The internal Seed Filter configuration, contained in the seed "
          "finder config, is not in internal units.");
    }

    using namespace UnitLiterals;
    SeedFinderConfig config = *this;
    config.isInInternalUnits = true;
    config.minPt /= 1_MeV;
    config.deltaRMin /= 1_mm;
    config.deltaRMax /= 1_mm;
    config.binSizeR /= 1_mm;
    config.deltaRMinTopSP /= 1_mm;
    config.deltaRMaxTopSP /= 1_mm;
    config.deltaRMinBottomSP /= 1_mm;
    config.deltaRMaxBottomSP /= 1_mm;
    config.deltaRMiddleMinSPRange /= 1_mm;
    config.deltaRMiddleMaxSPRange /= 1_mm;
    config.impactMax /= 1_mm;
    config.maxPtScattering /= 1_MeV;  // correct?
    config.collisionRegionMin /= 1_mm;
    config.collisionRegionMax /= 1_mm;
    config.zMin /= 1_mm;
    config.zMax /= 1_mm;
    config.rMax /= 1_mm;
    config.rMin /= 1_mm;
    config.deltaZMax /= 1_mm;

    config.zAlign /= 1_mm;
    config.rAlign /= 1_mm;

    config.toleranceParam /= 1_mm;

    return config;
  }
  SeedFinderConfig calculateDerivedQuantities() const {
    SeedFinderConfig config = *this;
    // calculation of scattering using the highland formula
    // convert pT to p once theta angle is known
    config.highland = 13.6 * std::sqrt(radLengthPerSeed) *
                      (1 + 0.038 * std::log(radLengthPerSeed));
    const float maxScatteringAngle = config.highland / minPt;
    config.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;
    return config;
  }
};

struct SeedFinderOptions {
  // location of beam in x,y plane.
  // used as offset for Space Points
  Vector2 beamPos{0 * UnitConstants::mm, 0 * UnitConstants::mm};
  // field induction
  float bFieldInZ = 2.08 * UnitConstants::T;

  // derived quantities
  float pTPerHelixRadius = std::numeric_limits<float>::quiet_NaN();
  float minHelixDiameter2 = std::numeric_limits<float>::quiet_NaN();
  float pT2perRadius = std::numeric_limits<float>::quiet_NaN();
  float sigmapT2perRadius = std::numeric_limits<float>::quiet_NaN();
  float multipleScattering2 = std::numeric_limits<float>::quiet_NaN();

  bool isInInternalUnits = false;

  SeedFinderOptions toInternalUnits() const {
    if (isInInternalUnits) {
      throw std::runtime_error(
          "Repeated conversion to internal units for SeedFinderOptions");
    }
    using namespace UnitLiterals;
    SeedFinderOptions options = *this;
    options.isInInternalUnits = true;
    options.beamPos[0] /= 1_mm;
    options.beamPos[1] /= 1_mm;

    options.bFieldInZ /= 1000. * 1_T;
    return options;
  }

  template <typename Config>
  SeedFinderOptions calculateDerivedQuantities(const Config& config) const {
    using namespace UnitLiterals;

    if (!isInInternalUnits) {
      throw std::runtime_error(
          "Derived quantities in SeedFinderOptions can only be calculated from "
          "Acts internal units");
    }
    SeedFinderOptions options = *this;
    // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
    // millimeter
    // TODO: change using ACTS units
    options.pTPerHelixRadius = 1_T * 1e6 * options.bFieldInZ;
    options.minHelixDiameter2 =
        std::pow(config.minPt * 2 / options.pTPerHelixRadius, 2) *
        config.helixCutTolerance;
    options.pT2perRadius =
        std::pow(config.highland / options.pTPerHelixRadius, 2);
    options.sigmapT2perRadius =
        options.pT2perRadius * std::pow(2 * config.sigmaScattering, 2);
    options.multipleScattering2 =
        config.maxScatteringAngle2 * std::pow(config.sigmaScattering, 2);
    return options;
  }
};

template <typename external_spacepoint_t, typename grid_t,
          typename platform_t = void*>
class SeedFinder {
 public:
  using Config = SeedFinderConfig<external_spacepoint_t>;

  using Options = SeedFinderOptions;

  struct State {
    // bottom space point
    std::vector<const external_spacepoint_t*> compatBottomSP{};
    std::vector<const external_spacepoint_t*> compatTopSP{};
    // contains parameters required to calculate circle with linear equation
    // ...for bottom-middle
    std::vector<LinCircle> linCircleBottom{};
    // ...for middle-top
    std::vector<LinCircle> linCircleTop{};

    // create vectors here to avoid reallocation in each loop
    std::vector<const external_spacepoint_t*> topSpVec{};
    std::vector<float> curvatures{};
    std::vector<float> impactParameters{};

    // managing seed candidates for SpM
    CandidatesForMiddleSp<const external_spacepoint_t> candidates_collector{};

    // managing doublet candidates
    boost::container::small_vector<Neighbour<grid_t>,
                                   detail::ipow(3, grid_t::DIM)>
        bottomNeighbours{};
    boost::container::small_vector<Neighbour<grid_t>,
                                   detail::ipow(3, grid_t::DIM)>
        topNeighbours{};

    // Mutable variables for Space points used in the seeding
    SpacePointMutableData spacePointMutableData{};
  };

  /// The only constructor. Requires a config object.
  /// @param config the configuration for the SeedFinder
  /// @param logger the ACTS logger
  explicit SeedFinder(const SeedFinderConfig<external_spacepoint_t>& config,
                      std::unique_ptr<const Logger> logger =
                          getDefaultLogger("Finder", Logging::Level::INFO));
  SeedFinder(SeedFinder<external_spacepoint_t, grid_t, platform_t>&&) noexcept =
      default;
  SeedFinder& operator=(SeedFinder<external_spacepoint_t, grid_t,
                                   platform_t>&&) noexcept = default;
  ~SeedFinder() = default;
  /**    @name Disallow default instantiation, copy, assignment */
  //@{
  SeedFinder() = default;
  SeedFinder(const SeedFinder<external_spacepoint_t, grid_t, platform_t>&) =
      delete;
  SeedFinder<external_spacepoint_t, grid_t, platform_t>& operator=(
      const SeedFinder<external_spacepoint_t, grid_t, platform_t>&) = delete;
  //@}

  /// Create all seeds from the space points in the three iterators.
  /// Can be used to parallelize the seed creation
  /// @param options frequently changing configuration (like beam position)
  /// @param state State object that holds memory used
  /// @param grid The grid with space points
  /// @param outputCollection Output container for the seeds in the group
  /// @param bottomSPs group of space points to be used as innermost SP in a seed.
  /// @param middleSPs group of space points to be used as middle SP in a seed.
  /// @param topSPs group of space points to be used as outermost SP in a seed.
  /// @param rMiddleSPRange range object containing the minimum and maximum r for middle SP for a certain z bin.
  /// @note Ranges must return pointers.
  /// @note Ranges must be separate objects for each parallel call.
  template <typename container_t, GridBinCollection sp_range_t>
    requires CollectionStoresSeedsTo<container_t, external_spacepoint_t, 3ul>
  void createSeedsForGroup(const SeedFinderOptions& options, State& state,
                           const grid_t& grid, container_t& outputCollection,
                           const sp_range_t& bottomSPs,
                           const std::size_t middleSPs,
                           const sp_range_t& topSPs,
                           const Range1D<float>& rMiddleSPRange) const;

 private:
  /// Given a middle space point candidate, get the proper radius validity range
  /// In case the radius range changes according to the z-bin we need to
  /// retrieve the proper range. We can do this computation only once, since
  /// all the middle space point candidates belong to the same z-bin
  /// @param spM space point candidate to be used as middle SP in a seed
  /// @param rMiddleSPRange range object containing the minimum and maximum r for middle SP for a certain z bin.
  std::pair<float, float> retrieveRadiusRangeForMiddle(
      const external_spacepoint_t& spM,
      const Range1D<float>& rMiddleSPRange) const;

  /// Iterates over dublets and tests the compatibility between them by applying
  /// a series of cuts that can be tested with only two SPs
  /// @param options frequently changing configuration (like beam position)
  /// @param grid spacepoint grid
  /// @param mutableData Container for mutable variables used in the seeding
  /// @param otherSPsNeighbours inner or outer space points to be used in the dublet
  /// @param mediumSP space point candidate to be used as middle SP in a seed
  /// @param linCircleVec vector containing inner or outer SP parameters after reference frame transformation to the u-v space
  /// @param outVec Output object containing top or bottom SPs that are compatible with a certain middle SPs
  /// @param deltaRMinSP minimum allowed r-distance between dublet components
  /// @param deltaRMaxSP maximum allowed r-distance between dublet components
  /// @param uIP minus one over radius of middle SP
  /// @param uIP2 square of uIP
  /// @param cosPhiM ratio between middle SP x position and radius
  /// @param sinPhiM ratio between middle SP y position and radius
  template <SpacePointCandidateType candidateType, typename out_range_t>
  void getCompatibleDoublets(
      const SeedFinderOptions& options, const grid_t& grid,
      SpacePointMutableData& mutableData,
      boost::container::small_vector<
          Neighbour<grid_t>, detail::ipow(3, grid_t::DIM)>& otherSPsNeighbours,
      const external_spacepoint_t& mediumSP,
      std::vector<LinCircle>& linCircleVec, out_range_t& outVec,
      const float deltaRMinSP, const float deltaRMaxSP, const float uIP,
      const float uIP2, const float cosPhiM, const float sinPhiM) const;

  /// Iterates over the seed candidates tests the compatibility between three
  /// SPs and calls for the seed confirmation
  /// @param SpM space point candidate to be used as middle SP in a seed
  /// @param options frequently changing configuration (like beam position)
  /// @param seedFilterState State object that holds memory used in SeedFilter
  /// @param state State object that holds memory used
  template <DetectorMeasurementInfo detailedMeasurement>
  void filterCandidates(const external_spacepoint_t& SpM,
                        const SeedFinderOptions& options,
                        SeedFilterState& seedFilterState, State& state) const;

 private:
  const Logger& logger() const { return *m_logger; }

  SeedFinderConfig<external_spacepoint_t> m_config;
  std::unique_ptr<const Logger> m_logger{nullptr};
};

}  // namespace Acts

#ifndef DOXYGEN
#include "Acts/Seeding/SeedFinder.ipp"
#endif

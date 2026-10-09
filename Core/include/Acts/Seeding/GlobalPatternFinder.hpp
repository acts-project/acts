// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/detail/GlobalPatternFinderAuxiliaries.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "Acts/Utilities/KDTree.hpp"

/** The Global Pattern Finder (GPF) is a fast reconstruction algorithm whose primary purpose is to 
 *  identify and collect detector hits belonging to the same particle trajectory. In contrast 
 *  to conventional tracking algorithms, the GPF does not attempt to estimate track parameters; 
 *  instead, it focuses on efficiently identifying compatible hit patterns that can subsequently 
 *  be used by downstream reconstruction algorithms.
 *
 *  The basic assumption of the GPF is that, to first order, the trajectory can be described by 
 *  a bending and a non-bending coordinate, where the bending coordinate describes the direction 
 *  in which the trajectory changes due to bending, while the non-bending coordinate remains 
 *  approximately constant over the relevant detector region. As a consequence, the detector can 
 *  be sectorized in the non-bending coordinate. This sectorization reduces the search space and 
 *  allows the GPF to efficiently identify compatible hits across the detector layers.
 *
 *  The GPF also abstracts the detector topology through a pattern topology, which describes how
 *  detector layers are organised into groups for pattern finding. The topology defines the
 *  ordering of layers, the grouping of layers, and the association of hits with the corresponding
 *  layer/group. This allows the pattern-finding algorithm to operate independently of the specific
 *  layer and grouping structure of the detector.
 *
 *  The precise definition of the bending and non-bending coordinates, as well as the corresponding 
 *  sectorization, is detector-dependent. The pattern-finding strategy itself is independent of 
 *  these specific coordinate choices.
 *  **** This first implementation is developed for a toroidal magnetic field where the bending
 *  coordinate is the polar angle theta, while the non-bending one is the azimuthal angle phi. ****
 *
 *  The definition of the concepts for the hit type, pattern topology, and sectorization, together
 *  with the main @ref PatternStateAux object used by the GPF, are provided in the GlobalPatternFinderAuxiliaries.
 */
namespace Acts::Experimental::detail {

template<typename Selector_t, typename Hit_t>
concept PatternSeedSelector = 
    GlobPatFinderHit<Hit_t> &&
    requires(const Selector_t& selector, 
             const Hit_t& hit) {
    { selector.goodForSeeding(hit) } -> std::same_as<bool>;
    { selector.thetaSearchWindow(hit) } -> std::same_as<double>;
};

template<typename Provider_t, typename Hit_t, typename Topology_t, typename Pattern_t>
concept OnlyPhiHitsProvider = 
    GlobPatFinderHit<Hit_t> &&
    PatternTopology<Topology_t, Hit_t> &&
    requires(const Provider_t& provider,
             const Pattern_t& pattern,
             const GeometryContext& gctx) {
    { provider.getPhiOnlyHits(pattern, gctx) } 
        -> std::same_as<std::array<std::vector<Hit_t>, Topology_t::nGroups>>;
};

template<GlobPatFinderHit Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
class GlobalPatternFinder {
  public:
    using PatternStateAux = detail::PatternStateAux<Hit_t, Sector_t, Topology_t>;
    using BeamspotInfo = PatternStateAux::BeamspotInfo;
    /** @brief Configuration object for the patter finder */     
    struct Config : PatternStateAux::Config {
        /********* Pattern bulding acceptance **********/ 
        /** @brief Number of standard deviations to consider for residual acceptance */
        double nResidualSigma {3.0};
        /** @brief Number of standard deviations to consider for phi acceptance */
        double nPhiSigma {3.0};
        /** @brief Size of theta window in radians to search for comapatible hits with a pattern, tailored to the target pt cutoff */
        double maxThetaSizeOverlap {0.05};
        /** @brief Residual uncertainty to consider the hit as low confidence */
        double lowConfidenceResSigma {50};
        /** @brief How much loosen the acceptance window when using beamspot or searching in a new group of the detector */
        double loosenAcceptanceFactor {2.};
        /********* Pile-up & Fake rate suppression *****/
        /** @brief Minimum number of strip layers in the bending direction required to accept a pattern */
        unsigned int minStripEtaLayers {3};
        /** @brief Minimum number of precision layers in the bending direction required to accept a pattern */
        unsigned int minPrecisionLayers {0};
        /** @brief Minimum number of phi layers required to accept a pattern */
        unsigned int minPhiLayers {1};
        /** @brief Minimum number of layer groups required to accept a pattern */
        unsigned int minGroups {2};
        /** @brief Minimum number of layers in a group to be considered a good group */
        unsigned int minGroupLayers {4};
        /** @brief Quality cut on pattern'mean squared normalized residual */
        double meanNormRes2Cut {0.2};
        /********* Pattern recovery power **************/
        /** @brief Maximum number of attempts to build a pattern from hits already used in existing patterns */
        unsigned int maxSeedAttempts {2};
        /** @brief Maximum number of missed candidate hits in different measurement layers in a group */
        unsigned int maxMissLayersInGroup {2};
        /********* Only-phi hit search *********/
        bool searchOnlyPhiHits {true};
        
    };
    /** @brief Abbreviation of the seed coordinates */
    enum class HitCoords : std::uint8_t {
        /** Sector coordinate of the associated spectrometer sector */
        eSector=0,
        /** Global Theta */
        eTheta=1,
        /** Number of coodinates */
        eNCoords=2
    };
        
    /** @brief Definition of the search tree class */
    using SearchTree_t = KDTree<toUnderlying(HitCoords::eNCoords), const Hit_t*, double, std::array, 50>;

    /** @brief Structure to hold the pattern result */
    struct OutputPattern {
        explicit OutputPattern(typename Sector_t::Index_t sector) 
            : sector{sector} {}
        /** @brief Vector of hits in the pattern */
        std::array<std::vector<const Hit_t*>, Topology_t::nGroups> hitsPerGroup{};
        /** @brief Vector of phi-only hits */
        std::vector<Hit_t> phiOnlyHits{};
        /** @brief Mean over eta hits of the square of their residual divided by residual uncertainty */
        double meanNormResidual2{0.};
        /** @brief Pattern phi, which is the phi of the bending plane where the pattern lies */
        double patPhi{0.};
        /** @brief Pattern theta, which is the value of the seed hit */
        double patTheta{0.};
        /** @brief Sector */
        Sector_t sector{static_cast<typename Sector_t::Index_t>(0u)};
        /** @brief Counts of precision / non-precision / phi layers  */
        uint8_t nPrecisionLayers{0u};
        uint8_t nTriggerLayers{0u};
        uint8_t nPhiLayers{0u};
    };

    /** @brief Standard constructor
     *  @param name: Name to be printed in the messaging
     *  @param config: Configuration parameters */
    explicit GlobalPatternFinder(Config&& config,
                                 std::unique_ptr<const Logger> logger = getDefaultLogger(
                                     "GlobalPatternFinder", Logging::Level::INFO));

    /** @brief Main methods steering the pattern finding. Given the space-point containers, it creates the search tree,  
     *         builds patterns in eta, attach compatible only-phi measurements, and convert PatternStateAuxs into GlobalPatterns
     *  @param gctx: Geometry context
     *  @param treeData: Data for the search tree
     *  @param seedSelector: Selector for the seed
     *  @param onlyPhiProvider: Provider for only-phi hits
     *  @param beamspotInfo: Beamspot information
     *  @return: Vector of found patterns */
    template<PatternSeedSelector<Hit_t> SeedSelector_t,
             OnlyPhiHitsProvider<Hit_t, Topology_t, PatternStateAux> OnlyPhiProvider_t>
    std::vector<OutputPattern> 
    findPatterns(const GeometryContext& gctx,
                 const SearchTree_t& treeData,
                 const SeedSelector_t& seedSelector,
                 const OnlyPhiProvider_t& onlyPhiProvider,
                 const BeamspotInfo& beamspotInfo) const;

  private:
    using PatternStateAuxVec = std::vector<PatternStateAux>;
    using OrderedHit = typename PatternStateAux::OrderedHit;
    using LineTestRes = typename PatternStateAux::LineTestRes;

    /** @brief Method steering the global pattern building in the bending plane.
     *  @param gctx: Geometry context
     *  @param orderedSpacepoints: Search tree with spacepoints ordered by their corresponding coordinates
     *  @param seedSelector: The seed selector to use for selecting the initial seed
     *  @param beamspotInfo: Information about the beam spot
     *  @return: resulting vector of PatternStateAuxs successfully built */
    template<PatternSeedSelector<Hit_t> SeedSelector_t>
    PatternStateAuxVec
    findPatternsInEta(const GeometryContext& gctx,
                      const SearchTree_t& orderedSpacepoints,
                      const SeedSelector_t& seedSelector,
                      const BeamspotInfo& beamspotInfo) const;
    /** @brief Main function controlling the development of patterns, including pattern branching when necessary. 
     *         It tests pattern compatibility of a set of active patterns (patterns produced from the same seed hit) against one 
     *         test hit. At the end, activePatterns contains the surviving patterns.
     *  @param startPatterns: Vector of active patterns to be extended
     *  @param endPatterns: Vector to store the surviving patterns after testing against the test hit.
     *  @param testHit: Hit to be tested against the patterns
     *  @param beamSpotInfo: Information about the beam spot, needed when the pattern line cannot be reliably defined from the pattern hits */
    void extendPatterns(const GeometryContext& gctx,
                        PatternStateAuxVec& startPatterns,
                        PatternStateAuxVec& endPatterns,
                        const OrderedHit& testHit,
                        const BeamspotInfo& beamSpotInfo) const;
    /** @brief Method checking line compatibility of a test hit against the pattern
     *  @param testHit: test hit information
     *  @param beamSpot: Beam spot position, needed to update the pattern line
     *  @return: result of the test, including the computed line residual and acceptance window */
    LineTestRes checkLineCompatibility(const GeometryContext& gctx,
                                       PatternStateAux& pat,
                                       const OrderedHit& testHit,
                                       const BeamspotInfo& beamSpot) const;
    /** @brief Method to check the phi compatibility of a test hit with the pattern
     *  @param gctx: geometry context
     *  @param pat: pattern to be checked
     *  @param hit: hit to be checked
     *  @return: true if the test hit is phi compatible with the pattern */
    bool isPhiCompatible(const GeometryContext& gctx,
                         const PatternStateAux& pat,
                         const Hit_t& hit) const;
    /** @brief Method to check if a pattern passes the quality cuts
     *  @param pattern: Pattern to be checked
     *  @return: true if the pattern passes the cuts, false otherwise */
    bool passPatternCuts(const PatternStateAux& pat) const;
    /** @brief Method to compare two patterns and define which one is better.
     *  @param a: first pattern
     *  @param b: second pattern
     *  @return: true if pattern a is better than pattern b, false otherwise */
    static bool isBetter(const PatternStateAux& a, 
                         const PatternStateAux& b);
    /** @brief Method to count the number of good groups in a pattern
     *  @param pat: pattern to be checked
     *  @return: number of good groups in the pattern */
    typename Topology_t::GroupIdx countGoodGroups(const PatternStateAux& pat) const;
    /** @brief Method to remove overlapping patterns
     *  @param toResolve: pattern to be resolved
     *  @return: resolved patterns */
    PatternStateAuxVec
    resolveOverlaps(PatternStateAuxVec& toResolve) const;
    /** @brief Method to add phi-only measurements to existing PatternStateAuxs
     *  @param gctx: Geometry context
     *  @param patterns: Vector of pattern states to which to add phi-only hits
     *  @return: Vector of added phi-only hits */
    template<OnlyPhiHitsProvider<Hit_t, Topology_t, PatternStateAux> OnlyPhiProvider_t>
    void addPhiOnlyHits(const GeometryContext& gctx,
                        const BeamspotInfo& beamSpot,
                        const OnlyPhiProvider_t& onlyPhiProvider,
                        PatternStateAuxVec& patterns) const;
    /** @brief Method to convert a PatternStateAux into a GlobalPattern object
     *  @param candidate: PatternStateAux to be converted
     *  @return: Converted GlobalPattern */
    OutputPattern
    convertToPattern(PatternStateAux&& candidate) const;
    /** @brief Method to convert a vector of PatternStateAuxs into GlobalPattern objects
     *  @param candidates: PatternStateAuxs to be converted
     *  @return: Vector of converted GlobalPatterns */
    std::vector<OutputPattern> 
    convertToPattern(PatternStateAuxVec&& candidates) const;

    /** @brief Global Pattern Recognition configuration */
    Config m_cfg;
    /** @brief Cut to the residual when we have 2 dof */
    double m_cut2Dof{std::sqrt(-2. * std::log(std::erfc(m_cfg.nResidualSigma / std::numbers::sqrt2)))};
    /** @brief Logger for the Global Pattern Finder */
    std::unique_ptr<const Acts::Logger> m_logger{};
    /// Reference to the logger object
    const Logger& logger() const { return *m_logger; }
};
}

#include "Acts/Seeding/GlobalPatternFinder.ipp"

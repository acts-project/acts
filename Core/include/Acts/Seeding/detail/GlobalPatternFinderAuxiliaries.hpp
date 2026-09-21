// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/CompositeSpacePoint.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts::Experimental::detail {

template <typename Hit_t>
concept GlobPatFinderHit = requires(const Hit_t hit,
                                    const GeometryContext& gctx,
                                    const Vector3& contractionVector,
                                    const Hit_t& otherHit) {

    /// Unit vector pointing to the next strip/straw in the plane or in case of a
    /// combined measurement, the complementary strip direction
    { hit.spacePoint() };
    requires CompositeSpacePointPtr<decltype(hit.spacePoint())>;
        
    /// Global position of the space point measurement. It's either the position of
    /// the wire or the position of the fired strip in global cooordinates. This method
    /// can either be implemented to return the vector by value or by reference.
    { hit.globalPosition(gctx) } -> std::convertible_to<const Vector3>;
    /// Orientation of the sensor, which is either the wire orientation or the
    /// strip orientation. Distortions along the sensor direction do not alter the
    /// track residual
    { hit.globalSensorDirection(gctx) } -> std::convertible_to<const Vector3>;
  
    /// Contraction of the hit position covariance along the given contraction vector.
    { hit.intrinsicVariance(gctx, contractionVector) } -> std::same_as<double>;
    /// Radius of the variance in [rad] of the space point measurement.
    { hit.phiVariance(gctx) } -> std::same_as<double>;
    /// Whether the hit is a precision measurement
    { hit.isPrecision() } -> std::same_as<bool>;
    /// Compare two hits for equality
    { hit == otherHit } -> std::same_as<bool>;
};

template<typename T, typename Hit_t>
concept PatternTopology = 
    GlobPatFinderHit<Hit_t> && 
    requires(const Hit_t& hit1, 
             const Hit_t& hit2,
             const Hit_t& hit) {
    typename T::GroupIdx;
    typename T::LayerIdx;
    requires std::unsigned_integral<typename T::GroupIdx>;
    requires std::unsigned_integral<typename T::LayerIdx>;

    { T::nGroups } -> std::convertible_to<typename T::GroupIdx>;
    { T::layerSorter(hit1, hit2) } -> std::same_as<bool>;
    { T::sameLayer(hit1, hit2) } -> std::same_as<bool>;
    { T::groupIndex(hit) } -> std::convertible_to<typename T::GroupIdx>;
};

template<typename Sector_t>
concept SectorType = requires(const Sector_t& sector,
                              const Sector_t& other,
                              const double phi) {
    typename Sector_t::Index_t;
    requires std::constructible_from<Sector_t, 
                                     typename Sector_t::Index_t>;

    { sector.sector() } -> std::same_as<typename Sector_t::Index_t>;
    { sector.phi() } -> std::same_as<double>;
    { sector.normalDir() } -> std::same_as<Vector3>;
    { sector.sectorSize() } -> std::same_as<double>;
    { sector.isNeighbour(other) } -> std::same_as<bool>;
    { sector.insideSector(phi) } -> std::same_as<bool>;
};

/** @brief Pattern state object storing pattern information during construction */
template<GlobPatFinderHit Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
struct PatternState {
    /// Configuration object of the residual calculator
    struct Config {
        /** @brief Number of standard deviations to consider for residual acceptance */
        double nResidualSigma {3.0};
        /** @brief Number of standard deviations to consider for phi acceptance */
        double nPhiSigma {3.0};
        /** @brief Minimum number of layers in a group to be considered a good group */
        unsigned int minGroupLayers {4};
        /********* Numerical stability *****************/
        /** @brief Minimum distance (in mm) between two hits for being used to 
         *         compute a reliable pattern line. Use the beamspot otherwise. */
        double minHitDistance4Line {40};
    };
    /** @brief Beamspot information */
    struct BeamspotInfo {
        /** @brief Beamspot position */
        Vector3 position{Vector3::Zero()};
        /** @brief Beamspot radius */
        double radius{0.};
        /** @brief Beamspot length */
        double length{0.};
    };
    /** @brief Small wrapper for candidate hits used to build patterns. This is needed
     *         because the global layer number cannot be defined globally, but it can be
     *         computed given a set of hits. */
    struct OrderedHit {
        /** @brief Pointer to the underlying hit */
        const Hit_t* hitPtr{nullptr};
        /** @brief Global measurement layer number */
        typename Topology_t::LayerIdx globLayer{0u};
        // Forward commonly used accessors for convenience
        const Hit_t* operator->() const { return hitPtr; }
        const Hit_t& operator*() const { return *hitPtr; }
        bool operator==(const OrderedHit& other) const { return *hitPtr == *other.hitPtr; }
        bool operator==(const Hit_t& other) const { return *hitPtr == other; }
        // Print and stream operator
        friend std::ostream& operator<<(std::ostream& ostr, const OrderedHit& c) {
            c.print(ostr);
            return ostr;
        }
        void print(std::ostream& ostr) const;
    };
    /** @brief: Enum for possible outcomes of pattern line compatibility test */        
    enum class LineTestDecision : std::int8_t{
        /** @brief Test successfull, add hit to pattern */
        eAddHit,
        /** @brief Test successfull with multiple pattern hits on same layer, branch the pattern */
        eBranchPattern,
        /** @brief Test failed, discard the hit */
        eRejectHit
    };
    /** @brief : Small struct to encapsulate the result of the line compatibility test */
    struct LineTestRes  {
        double residual{0.};
        double sigma{0.};
        LineTestDecision result {LineTestDecision::eRejectHit};
    };
    /** @brief Constructor taking the seed information 
        *  @param seed: seed hit
        *  @param expSector: sector index of the seed hit
        *  @param cfg: pointer to configuration object
        *  @param logger: pointer to messaging object */
    explicit PatternState(const GeometryContext& gctx,
                          const OrderedHit& seed,
                          const typename Sector_t::Index_t expSector,
                          const Config* cfg,
                          const Logger* logger);
    /** @brief Move constructor
        *  @param other: other pattern state to move from */
    PatternState(PatternState&& other) noexcept = default;
    /** @brief Move assignment operator
        *  @param other: other pattern state to move from */
    PatternState& operator=(PatternState&& other) noexcept = default;
    /** @brief Copy constructor
        *  @param other: other pattern state to copy from */
    PatternState(const PatternState& other) = default;
    /** @brief Copy assignment operator
        *  @param other: other pattern state to copy from */
    PatternState& operator=(const PatternState& other) = default;
    /** @brief Default destructor */
    ~PatternState() =default;

    /** @brief Add a hit to the pattern and update the internal state
        *  @param hit: hit to be added
        *  @param residual: residual of the hit
        *  @param resSigma: residual uncertainty of the hit */
    void addHit(const GeometryContext& gctx,
                const OrderedHit& hit,
                const double residual,
                const double resSigma);
    /** @brief Overwrite the hits on the last layer with the new one
        *  @param newHit: new hit to replace with
        *  @param newResidual: residual of the new hit 
        *  @param newResSigma: residual uncertainty of the new hit
        *  @param beamSpot: needed to update line parameters */
    void overWriteHit(const GeometryContext& gctx,
                      const OrderedHit& newHit,
                      const double newResidual,
                      const double newResSigma);
    /** @brief Method to compute the residual of a test hit against the pattern line
        *  @param testHit: test hit information
        *  @param beamSpot: position of the beam spot
        *  @return: Test result holding the residual and acceptance window. The decision is set later. */
    LineTestRes computeLineResidual(const GeometryContext& gctx,
                                    const OrderedHit& testHit,
                                    const BeamspotInfo& beamSpot) const;
    /** @brief Project a certain hit position onto the bending plane where the pattern is defined. 
        *         The hit is moved along the sensor direction if it does not measure phi, 
        *         or is rotated around the Z axis if it does.
        *  @param hit: hit whose position is to be projected
        *  @return: projected position */
    Vector3 projToPhiPlane(const GeometryContext& gctx,
                           const Hit_t& hit) const;
    /** @brief Method to check the phi compatibility of a test hit with a given pattern
        *  @param gctx: geometry context
        *  @param hit: hit to be checked
        *  @return: true if the test hit is phi compatible with the pattern, false otherwise */
    bool isPhiCompatible(const GeometryContext& gctx, 
                         const Hit_t& hit) const;
    /** @brief Check wheter a hit is present in the pattern
        *  @param hit: hit to be checked
        *  @return: boolean indicating if the hit is in the pattern */
    bool isInPattern(const Hit_t& hit) const;
    /** @brief Move the line anchor hit given a reference hit. The anchor is defined
                as the closest hit in the closest group to the referece hit, */
    void moveLineAnchorHit(const OrderedHit& refHit);
    /** @brief Update the line parameters based on the current hits
        *  @param beamSpot: position of the beam spot, needed when there are not enough hits */
    void updateLineParameters(const GeometryContext& gctx,
                              const BeamspotInfo& beamSpot);
    /** @brief Helper method to update the pattern phi and bending plane normal */
    void updatePatternPhi(const GeometryContext& gctx);
    /** @brief Return the mean normalized residual squared */
    double getMeanResidual2() const;
    /** @brief Method returning the number of groups
     *  @param onlyGoodGroups: flag to indicate if only good groups should be counted,
     *         i.e. having a minimum number of hits */
    typename Topology_t::GroupIdx countGroups(const bool onlyGoodGroups) const;
    /** @brief Return the number of layers in bending coordinate */
    typename Topology_t::LayerIdx nBendingLayers() const;
    
    /** @brief Return the logger */
    const Logger& logger() const {
        return *m_logger;
    }

    /** @brief Pointer to cfg option */
    const Config* cfg{nullptr};
    /** @brief Logger */
    const Logger* m_logger{nullptr};
    /** @brief Last inserted hit. Needed to speed-up lookup */
    OrderedHit lastInsertedHit{};
    /** @brief Last hit in the second-to-last layer */
    OrderedHit prevLayerHit{};
    /** @brief Line anchor hit */
    OrderedHit lineAnchorHit{};
    /** @brief Normal vector to the bending plane where the pattern lies */
    Vector3 bendPlaneNorm{Vector3::Zero()};
    /** @brief Position and direction of the pattern line. Both are constructed to be
        *         within the bending plane of the pattern. We store them as vectors to
        *         facilitate vector operations and avoid constructing them repeatedly. */
    Vector3 linePos{Vector3::Zero()};
    Vector3 lineDir{Vector3::Zero()};
    /** @brief Distance between the two points defining the pattern line */
    double leverArm{0.};
    /** @brief Mean over eta hits of the square of their residual divided by residual uncertainty */
    double meanNormResidual2{0.};
    /** @brief Residual & residual uncertainty of the last inserted hit (needed when replacing a hit) */
    double lastResidual{0.};
    double lastResSigma{0.};
    /** @brief Pattern phi, which is the phi of the bending plane where the pattern lies */
    double patPhi{0.};
    /** @brief Pattern theta, which is the value of the seed hit */
    double patTheta{0.};
    /** @brief Covariance of the pattern phi */
    double patPhiCov{0.};
    /** @brief Sector */
    Sector_t expSect{static_cast<typename Sector_t::Index_t>(0u)};
    /** @brief Counts of precision / non-precision / phi layers  */
    typename Topology_t::LayerIdx nPrecisionLayers{0u};
    typename Topology_t::LayerIdx nTriggerLayers{0u};
    typename Topology_t::LayerIdx nPhiLayers{0u};
    /** @brief Flag to indicate if the pattern has been finalized */
    bool isFinalized{false};
    /** @brief Flag to indicate if the pattern is overlapping with another one, used during overlap removal */
    bool isOverlap{false};
    /** @brief Whether we used the beamspot to compute the line parameters */
    bool useBeamspot{false};
    /** @brief Whether we need to update the pattern line the next time we find a hit in a new layer */
    bool needLineUpdate{false};
    
    /** @brief Map collection of hits per station. A pattern is determined by the hits belonging to it. */
    std::array<std::vector<OrderedHit>, Topology_t::nGroups> hitsPerGroup{};
    /** @brief Array holding phi-only hits */
    std::vector<Hit_t> phiOnlyHits{};

    /** @brief Patterns are considered identical if they have the same hit content. However, map comparison is very expensive */
    bool operator==(const PatternState& other) const = delete;
    /** @brief Print the pattern candidate */
    void print(std::ostream& ostr, bool detailed) const;
    /** @brief A view of the pattern state for printing purposes */
    struct PatternPrintView {
        /** @brief The pattern state to be printed */
        const PatternState& pat;
        /** @brief Whether to print detailed information */
        bool detailed = false;

        /** @brief Print the pattern print view */
        friend std::ostream& operator<<(std::ostream& os, const PatternPrintView& v) {
            v.pat.print(os, v.detailed);
            return os;
        }
    };
    /** @brief Print the pattern state with brief information */
    static PatternPrintView brief(const PatternState& p);
    /** @brief Print the pattern state with detailed information */
    static PatternPrintView detailed(const PatternState& p);
};

}
#include "Acts/Seeding/detail/GlobalPatternFinderAuxiliaries.ipp"
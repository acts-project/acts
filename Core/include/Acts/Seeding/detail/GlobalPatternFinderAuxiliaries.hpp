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
#include "Acts/Definitions/Units.hpp"

/** Auxiliary concepts and data structures used by the Global Pattern Finder.
 *  
 * This file defines the concepts used to abstract the detector-specific hit type,
 * pattern topology, and sectorization, together with the main PatternStateAux object
 * used internally by the GPF. */

namespace Acts::Experimental::detail {
/** @brief Definition of the measured directions for strip measurements. The first
 *  is the primary strip measurement direction and is by definition orthogonal to 
 *  the sensor direction. The second direction is the secondary strip measurement  
 *  direction and is provided only when it differs from the sensor direction. */
using StripMeasurementDirections =
    std::pair<Vector3, std::optional<Vector3>>;

/** @brief Concept definition for global pattern finder hits. This concept
 *  is defined as the representation of @ref CompositeSpacePoint in global
 *  coordinates, with additional methods needed for the global pattern finder.
 *  The transformation of represetation from in-station to global detector
 *  coordinates introduces a dependency on the geometry aligment, represented  
 *  by the @ref GeometryContext. */
template <typename Hit_t>
concept GlobPatFinderHit = requires(const Hit_t hit,
                                    const GeometryContext& gctx,
                                    const Vector3& contractionVector,
                                    const Hit_t& otherHit) {
    /// Retrieve the underlying composite space point. This is a proxy to all
    /// the properties of the composite space point.
    { hit.spacePoint() };
    requires CompositeSpacePointPtr<decltype(hit.spacePoint())>;
    /// Cached global position of the space point measurement. It's the position 
    /// of the space point in global coordinates.
    { hit.globalPosition(gctx) } -> std::convertible_to<Vector3>;
    /// Orientation of the sensor, which is either the wire or strip orientation. 
    /// When the spacepoint is built from two strip measurements, the sensor
    /// orientation refers to the measurement in the bending coordinate. This method   
    /// can either be implemented to return the vector by value or by reference.
    { hit.globalSensorDirection(gctx) } -> std::convertible_to<Vector3>;
    /// Return the measured directions for strip measurements. The first is orthogonal to 
    /// the sensor direction and is always present. The second is provided only when it 
    /// differs from the sensor direction.
    { hit.stripMeasurementDirections(gctx) } -> std::same_as<StripMeasurementDirections>;
    /// Whether the hit is a precision measurement
    { hit.isPrecision() } -> std::same_as<bool>;
    /// Whether two hits represent overlapping detector information. This relation may 
    /// be less restrictive than operator== and is used when determining whether 
    /// reconstructed patterns share detector information.
    { hit.overlaps(otherHit) } -> std::same_as<bool>;
};
/** @brief Equality operator for hits. */
template <GlobPatFinderHit Hit_t>
bool operator==(const Hit_t& lhs, const Hit_t& rhs) {
    return lhs.spacePoint() == rhs.spacePoint();
}

/** @brief Concept defining the detector topology used by the pattern finder. The
 *  detector is assumed to be organized in layers, which are grouped into a number 
 *  of groups. The topology describes how detector layers are ordered and organized 
 *  into groups, and how hits are associated with the corresponding layer/group. */
template<typename T, typename Hit_t>
concept PatternTopology = 
    GlobPatFinderHit<Hit_t> && 
    requires(const Hit_t& hit1,
             const Hit_t& hit2) {
    /// Type definitions for the group and layer index types. These types
    /// must be unsigned integral types. The layer index type is assumed 
    /// to be able to represent all the layers that can be crossed by a particle.
    typename T::GroupIdx;
    typename T::LayerIdx;
    requires std::unsigned_integral<typename T::GroupIdx>;
    requires std::unsigned_integral<typename T::LayerIdx>;
    
    /// Return the number of groups in the topology
    { T::nGroups } -> std::same_as<const typename T::GroupIdx&>;
    /// Sort two hits according by measurement layer and, if they are on 
    /// the same layer, by their position along the bending coordinate.
    { T::layerSorter(hit1, hit2) } -> std::same_as<bool>;
    /// Check if two hits are on the same measurement layer.
    { T::sameLayer(hit1, hit2) } -> std::same_as<bool>;
    /// Return the group index of a hit.
    { T::groupIndex(hit1) } -> std::same_as<typename T::GroupIdx>;
};
/** @brief Sectorization of the detector in the non-bending coordinate. 
 *  In this first implementation of GPF with toroidal azimuthal field,
 *  the sectors are defined along the phi direction. */
template<typename Sector_t>
concept SectorType = requires(const Sector_t& sector,
                              const Sector_t& other,
                              const double phi) {
    /// Type definition for the sector index type.
    typename Sector_t::Index_t;
    /// The sector type must be constructible from the sector index.
    requires std::constructible_from<Sector_t, 
                                     typename Sector_t::Index_t>;
    /// Returns the sector index of the sector
    { sector.sector() } -> std::same_as<typename Sector_t::Index_t>;
    /// Returns the phi coordinate of the sector center
    { sector.phi() } -> std::same_as<double>;
    /// Returns the angular size of the sector in [rad]
    { sector.sectorSize() } -> std::same_as<double>;
    /// Returns whether the sector is a neighbour of another sector.
    /// This helper incorporate the wrap-around of the sectorization.
    { sector.isNeighbour(other) } -> std::same_as<bool>;
    /// Returns whether a given phi coordinate is inside the sector.
    { sector.insideSector(phi) } -> std::same_as<bool>;
};
using namespace Acts::UnitLiterals;
/** @brief Pattern state object storing pattern information during construction */
template<GlobPatFinderHit Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
struct PatternStateAux {
    /** @brief Configuration object of the residual calculator */
    struct Config {
        /** @brief Minimum distance (in mm) between two hits for being used to 
         *         compute a reliable pattern line. Use the beamspot otherwise. */
        double minHitDistance4Line {40._mm};
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
    /** @brief Small wrapper for hits used to build patterns. In general, due to detector 
     *         complexity, the layer number cannot be defined globally, but it can be 
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
    enum class LineTestDecision : uint8_t{
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
        uint8_t nDof{1u};
        LineTestDecision result {LineTestDecision::eRejectHit};
    };
    /** @brief Enum for the bending plane parameter covariance */
    enum class BendPlaneCov : std::int8_t {
        ePhiPhi, // Variance of the phi parameter
        eSS,     // Variance of the displacement parameter
        ePhiS    // Covariance of the phi parameter with the displacement parameter
    };
    using BendPlaneCov_t = std::array<double, 3>;
    /** @brief Constructor taking the seed information 
     *  @param gctx: geometry context
     *  @param seed: seed hit
     *  @param sector: sector index of the seed hit
     *  @param cfg: pointer to configuration object
     *  @param logger: pointer to messaging object */
    explicit PatternStateAux(const GeometryContext& gctx,
                             const OrderedHit& seed,
                             const typename Sector_t::Index_t sector,
                             const BeamspotInfo& beamSpot,
                             const Config* cfg,
                             const Logger* logger);
    /** @brief Move constructor
     *  @param other: other pattern state to move from */
    PatternStateAux(PatternStateAux&& other) noexcept = default;
    /** @brief Move assignment operator
     *  @param other: other pattern state to move from */
    PatternStateAux& operator=(PatternStateAux&& other) noexcept = default;
    /** @brief Copy constructor
     *  @param other: other pattern state to copy from */
    PatternStateAux(const PatternStateAux& other) = default;
    /** @brief Copy assignment operator
     *  @param other: other pattern state to copy from */
    PatternStateAux& operator=(const PatternStateAux& other) = default;
    /** @brief Default destructor */
    ~PatternStateAux() =default;

    /** @brief Add a hit to the pattern and update the internal state
     *  @param gctx: geometry context
     *  @param hit: hit to be added
     *  @param residual: residual of the hit
     *  @param resSigma: residual uncertainty of the hit */
    void addHit(const GeometryContext& gctx,
                const OrderedHit& hit,
                const double residual,
                const double resSigma,
                const BeamspotInfo& beamSpot);
    /** @brief Overwrite the hits on the last layer with the new one
     *  @param gctx: geometry context
     *  @param newHit: new hit to replace with
     *  @param newResidual: residual of the new hit 
     *  @param newResSigma: residual uncertainty of the new hit */
    void overWriteHit(const GeometryContext& gctx,
                      const OrderedHit& newHit,
                      const double newResidual,
                      const double newResSigma,
                      const BeamspotInfo& beamSpot);
    /** @brief Compute the contribution of the intrinsic covariance of the hit to the residual
     *  @param gctx: geometry context
     *  @param hit: hit for which to compute the intrinsic variance
     *  @param contractionVector: vector used to contract the covariance
     *  @param isProjected: flag indicating if the hit is projected
     *  @return: the intrinsic variance of the hit */
    static double intrinsicVariance(const GeometryContext& gctx,
                                    const Hit_t& hit,
                                    const Vector3& contractionVector,
                                    bool isProjected);
    /** @brief Method to compute the residual of a test hit against the pattern line
     *  @param gctx: geometry context
     *  @param testHit: test hit information
     *  @param beamSpot: position of the beam spot
     *  @return: Test result holding the residual and acceptance window. The decision is 
     *           set later from the GPF. */
    LineTestRes computeLineResidual(const GeometryContext& gctx,
                                    const OrderedHit& testHit,
                                    const BeamspotInfo& beamSpot) const;
    /** @brief Project a certain hit position onto the bending plane where the pattern is defined. 
     *         The hit is moved along the sensor direction.
     *  @param hit: hit whose position is to be projected
     *  @return: projected position */
    Vector3 projToPhiPlane(const GeometryContext& gctx,
                           const Hit_t& hit) const;
    /** @brief Check whether a hit is present in the pattern
     *  @param hit: hit to be checked
     *  @return: true if the hit is in the pattern*/
    bool isInPattern(const Hit_t& hit) const;
    /** @brief Move the line anchor hit given a reference hit. The anchor is defined
     *         as the closest hit in the closest group to the referece hit, */
    void moveLineAnchorHit(const OrderedHit& refHit);
    /** @brief Update the line parameters based on the current hits
     *  @param beamSpot: position of the beam spot, needed when there are not enough hits */
    void updateLineParameters(const GeometryContext& gctx,
                              const BeamspotInfo& beamSpot);
    /** @brief Helper method to update the pattern phi and bending plane normal
     *  @param gctx: geometry context
     *  @param beamSpot: position of the beam spot */
    void updatePatternPhi(const GeometryContext& gctx,
                          const BeamspotInfo& beamSpot);
    /** @brief Signed angle under which the pattern sees a point, relative to its bending plane.
     *  @param pos: position to check
     *  @return: the angle between the position and the bending plane */
    double angleToBendPlane(const Vector3& pos) const;
    /** @brief Return the pattern's phi at a given radius */
    double phiAtRadius(const double R) const;
    /** @brief Return the mean normalized residual squared */
    double getMeanResidual2() const;
    /** @brief Method returning the number of groups
     *  @return: number of groups in the pattern */
    typename Topology_t::GroupIdx countGroups() const;
    /** @brief Return the number of layers in bending coordinate 
     *  @return: number of layers in the bending coordinate */
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
    /** @brief Signed offset of the bending plane from the origin */
    double patPhiOffset{0.};
    /** @brief Covariance of the bending plane */
    BendPlaneCov_t bendPlaneCov{};
    /** @brief Pattern theta, which is the value of the seed hit */
    double patTheta{0.};
    /** @brief Sector */
    Sector_t sector{static_cast<typename Sector_t::Index_t>(0u)};
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

    /** @brief Print the covariance matrix */
    static std::string print(const BendPlaneCov_t& cov) {
        std::ostringstream os;
        os << "[varPhi: " << cov[static_cast<std::size_t>(BendPlaneCov::ePhiPhi)]
           << ", varS: " << cov[static_cast<std::size_t>(BendPlaneCov::eSS)]
           << ", covPhiS: " << cov[static_cast<std::size_t>(BendPlaneCov::ePhiS)] << "]";
        return os.str();
    }

    /** @brief Patterns are considered identical if they have the same hit content. However, map comparison is very expensive */
    bool operator==(const PatternStateAux& other) const = delete;
    /** @brief Print the pattern candidate */
    void print(std::ostream& ostr, bool detailed) const;
    /** @brief A view of the pattern state for printing purposes */
    struct PatternPrintView {
        /** @brief The pattern state to be printed */
        const PatternStateAux& pat;
        /** @brief Whether to print detailed information */
        bool detailed = false;

        /** @brief Print the pattern print view */
        friend std::ostream& operator<<(std::ostream& os, const PatternPrintView& v) {
            v.pat.print(os, v.detailed);
            return os;
        }
    };
    /** @brief Print the pattern state with brief information */
    static PatternPrintView brief(const PatternStateAux& p);
    /** @brief Print the pattern state with detailed information */
    static PatternPrintView detailed(const PatternStateAux& p);
};

}
#include "Acts/Seeding/detail/GlobalPatternFinderAuxiliaries.ipp"
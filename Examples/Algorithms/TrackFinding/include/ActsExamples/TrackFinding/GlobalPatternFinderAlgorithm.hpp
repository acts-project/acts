// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Seeding/GlobalPatternFinder.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/OstreamFormatter.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <array>
#include <compare>
#include <cstdint>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

namespace ActsExamples {

/// @brief Sector of the muon spectrometer extended by the overlaps between
///        adjacent sectors. The spectrometer has a 16-fold symmetry in phi;
///        the expanded sector index is twice the spectrometer sector plus the
///        projector, which is -1 for the overlap with the previous sector, 0
///        for the sector centre and 1 for the overlap with the next sector.
class ExpandedSector {
 public:
  /// @brief Type of the index of the expanded sector
  using Index_t = std::int8_t;
  /// @brief Part of a spectrometer sector an expanded sector covers
  enum class SectorProjector : std::int8_t {
    /// Overlap with the previous sector
    leftOverlap = -1,
    /// Sector centre
    center = 0,
    /// Overlap with the next sector
    rightOverlap = 1
  };
  /// @brief Name of a projector
  /// @param proj Projector of interest
  /// @return The name of the enumerator
  static std::string toString(const SectorProjector proj);
  /// @brief Print a projector
  /// @param ostr Stream to print to
  /// @param proj Projector to print
  /// @return The stream
  friend std::ostream& operator<<(std::ostream& ostr,
                                  const SectorProjector proj) {
    return ostr << (toString(proj));
  }
  /// @brief Print an expanded sector
  /// @param ostr Stream to print to
  /// @param sec Expanded sector to print
  /// @return The stream
  friend std::ostream& operator<<(std::ostream& ostr,
                                  const ExpandedSector& sec) {
    return sec.toString(ostr);
  }
  /// @brief Constructor from a spectrometer sector and a projector
  /// @param msSector Number of the spectrometer sector [1-16]
  /// @param proj Part of the sector: the centre or the overlap with the
  ///             previous or the next sector
  explicit ExpandedSector(const unsigned msSector, const SectorProjector proj);
  /// @brief Constructor from a phi angle, which is assigned to the sector or
  ///        the overlap it lies in
  /// @param phi Angle in [-pi, pi]
  explicit ExpandedSector(const double phi);
  /// @brief Constructor from an expanded sector index
  /// @param expSector Expanded sector index
  explicit ExpandedSector(const Index_t expSector);
  /// @brief Define the ordering operator
  bool operator<(const ExpandedSector& other) const;
  /// @brief Define the equal operator
  bool operator==(const ExpandedSector& other) const;
  /// @brief Define the unequal operator
  bool operator!=(const ExpandedSector& other) const;
  /// @brief Returns the spectrometer sector of the expanded sector
  /// @note An overlap is stored as the right overlap of the lower sector,
  ///       hence the left overlap of a sector returns the previous sector
  unsigned msSector() const;
  /// @brief Returns the spectrometer sector on the other side of the overlap,
  ///        or msSector() for a sector centre
  unsigned adjacentMsSector() const;
  /// @brief Returns the projector in the corresponding MS sector
  SectorProjector projector() const;
  /// @brief Returns the expanded sector number
  Index_t sector() const;
  /// @brief Returns the phi angle of the expanded sector
  double phi() const;
  /// @brief Returns the vector that is normal to the plane spanned
  ///         by the expanded sector
  Acts::Vector3 normalDir() const;
  /// @brief Returns true if the expanded sector is a neighbour of the other
  bool isNeighbour(const ExpandedSector& other) const;
  /// @brief Check if a given phi is within the expanded sector
  bool insideSector(const double phi) const;
  /// @brief Return the expanded sector (half) size
  double sectorSize() const;

 private:
  /// @brief Print the expanded sector
  /// @param ostr Stream to print to
  /// @return The stream
  std::ostream& toString(std::ostream& ostr) const;
  /// @brief Expanded sector index
  std::int8_t m_sector{0};
};

static_assert(Acts::Experimental::detail::SectorType<ExpandedSector>);

/// @brief Station index of the muon spectrometer: the inner, middle, outer and
///        extended stations of the barrel (B) and of the endcap (E)
enum class StIndex : std::int8_t {
  StUnknown = -1,
  BI,
  BM,
  BO,
  BE,
  EI,
  EM,
  EO,
  EE,
  StIndexMax
};

/// @brief Station layer of the muon spectrometer, common to barrel and endcap
enum class LayerIndex : std::int8_t {
  LayerUnknown = -1,
  Inner,
  Middle,
  Outer,
  /// Endcap extended station (EE)
  Extended,
  /// Barrel extended station (BEE)
  BarrelExtended,
  LayerIndexMax
};

/// @brief Station index of a chamber
/// @param index Station name of the chamber
/// @return The station index, StUnknown for an undefined name
StIndex toStationIndex(MuonSpacePoint::MuonId::StationName index);

/// @brief Station layer of a station index
/// @param index Station index of interest
/// @return The station layer, LayerIndexMax for an unknown station
LayerIndex toLayerIndex(StIndex index);

/// @brief Returns true if the station index points to a barrel chamber
/// @param index Station index of interest
bool isBarrel(const StIndex index);

/// @brief Returns whether a space point is a precision hit: an MDT, a
///        micromegas or an sTGC strip measurement
/// @param hit Space point of interest
bool isPrecisionHit(const MuonSpacePoint& hit);

/// @brief Muon space point as a hit of the global pattern finder: the space
///        point together with its global position, station and measurement
///        surface
struct HitPayload {
  /// @brief Constructor
  /// @param gctx Geometry context
  /// @param sp Space point
  /// @param parentBucket Space point bucket holding the space point
  /// @param localToGlobal Transformation from the sector frame of the space
  ///                      point into the global frame
  /// @param measSurface Measurement surface of the space point
  explicit HitPayload(const Acts::GeometryContext& gctx,
                      const MuonSpacePoint* sp,
                      const MuonSpacePointBucket* parentBucket,
                      const Acts::Transform3& localToGlobal,
                      const Acts::Surface* measSurface);
  /// @brief Retrieve the space point
  const MuonSpacePoint& spacePoint() const { return *underlyingSp; }
  /// @brief Retrieve the global position
  const Acts::Vector3& globalPosition(const Acts::GeometryContext& gctx) const;
  /// @brief Retrieve the sensor direction
  Acts::Vector3 globalSensorDirection(const Acts::GeometryContext& gctx) const;
  /// @brief Directions of the bending and the non-bending measurement of a
  ///        strip hit. The non-bending direction is only given for tilted
  ///        strips; otherwise it is the sensor direction
  Acts::Experimental::detail::StripMeasurementDirections
  stripMeasurementDirections(const Acts::GeometryContext& gctx) const;
  /// @brief Returns whether the hit is a precision hit
  bool isPrecision() const { return isPrecisionFlag; }

  /// @brief Global position
  Acts::Vector3 position{Acts::Vector3::Zero()};
  /// @brief Pointer to the underlying hit
  const MuonSpacePoint* underlyingSp{nullptr};
  /// @brief Pointer to the parent bucket
  const MuonSpacePointBucket* bucket{nullptr};
  /// @brief Surface the measurement was recorded on. It defines the local axes
  ///         of the hit, i.e. the measuring, the sensor and the normal
  ///         direction
  const Acts::Surface* surface{nullptr};
  /// @brief Strip angle when the strips are non-orthogonal
  double stripAngle{0.};
  /// @brief Station index
  StIndex station{toStationIndex(spacePoint().id().msStation())};
  /// @brief Layer number in the sector frame
  std::uint8_t locLayer{
      static_cast<std::uint8_t>(spacePoint().id().detLayer())};
  /// @brief Is precision hit
  bool isPrecisionFlag{isPrecisionHit(spacePoint())};
  /// @brief Are the strips non-orthogonal
  bool nonOrthogonalStrips{false};
  /// @brief Equal operator: it compares the underlying hit
  bool operator==(const HitPayload& other) const;
  /// @brief Returns whether the two hits stem from the same measurement, which
  ///         is weaker than the equal operator: one measurement is held by
  ///         several space points
  bool overlaps(const HitPayload& other) const;
};
static_assert(Acts::Experimental::detail::GlobPatFinderHit<HitPayload>);

/// @brief Layer order and grouping of the muon spectrometer hits for the
///        global pattern finder. A group is a station index (BI, BM, ...);
///        the layers are ordered from the inside of the spectrometer outwards.
class PatternTopology {
 public:
  /// @brief Type of the index of the station layer
  using LayerIdx = std::uint16_t;
  /// @brief Type of the index of the group
  using GroupIdx = std::uint8_t;
  /// @brief Number of groups in a pattern
  static constexpr GroupIdx nGroups{
      static_cast<GroupIdx>(Acts::toUnderlying(StIndex::StIndexMax))};
  /// @brief Order of the layers of two hits, going from the inside of the
  ///        spectrometer outwards. Hits of one station are ordered by their
  ///        layer, hits of different stations by their station layer (inner,
  ///        extended, middle, outer). Two hits are on the same layer if they
  ///        are on the same layer of one chamber, or at the same radius
  ///        (barrel) or |z| (endcap) in the same station index
  /// @param hit1 First hit
  /// @param hit2 Second hit
  /// @return less if hit1 is on an inner layer, greater if on an outer layer,
  ///         equivalent if both hits are on the same layer
  /// @throws std::runtime_error for one hit in BO and the other in EO, which a
  ///         pattern cannot hold
  static std::weak_ordering compareLayers(const HitPayload& hit1,
                                          const HitPayload& hit2);
  /// @brief Group of a hit, its station index
  /// @param hit Hit of interest
  /// @return The group index of the hit
  static GroupIdx groupIndex(const HitPayload& hit);
};
static_assert(
    Acts::Experimental::detail::PatternTopology<PatternTopology, HitPayload>);

/// @brief Global pattern finder on the muon space points
using GlobalPatternFinder_t =
    Acts::Experimental::detail::GlobalPatternFinder<HitPayload, ExpandedSector,
                                                    PatternTopology>;

/// @brief Search tree over the hits of an event, in theta and expanded sector
using SearchTree_t = GlobalPatternFinder_t::SearchTree_t;

/// @brief Structure to hold the search tree data
struct SearchTreeData {
  /// @brief The hit payloads the tree points to
  std::vector<HitPayload> hitPayloads;
  /// @brief The search tree
  SearchTree_t tree;
};

/// @brief Chooses the hits the global pattern finder starts patterns from and
///        the theta windows it searches around each of them
class SeedSelector {
 public:
  /// @brief Configuration of the seed selection. The theta half-windows [rad]
  ///        are given per pair of seed station and hit station, for three bins
  ///        of the seed's |eta|; pairs left out are 0, so the hits of that
  ///        station are never candidates of the seed.
  struct Config {
    /// @brief Whether MDT hits may seed patterns
    bool seedFromMdt{false};
    /// @brief Whether hits in the inner station layer may seed patterns
    bool seedFromInner{false};
    /// @brief Upper |eta| edge of the barrel bin
    double etaBarrelMax{0.9};
    /// @brief Lower |eta| edge of the endcap bin; the transition bin lies
    ///        between
    double etaEndcapMin{1.45};

    /// @brief Barrel bin: two hits in the same station
    double barrelSameStation{0.9 * Acts::UnitConstants::degree};
    /// @brief Barrel bin: BI-BM and BM-BO
    double barrelAdjacent{1.29 * Acts::UnitConstants::degree};
    /// @brief Barrel bin: BI-BO
    double barrelSkip{1.861 * Acts::UnitConstants::degree};

    /// @brief Transition bin: two hits in the same station
    double transSameStation{0.888 * Acts::UnitConstants::degree};
    /// @brief Transition bin: BI-BM and BM-BO
    double transBarrelAdj{1.186 * Acts::UnitConstants::degree};
    /// @brief Transition bin: BI-BO
    double transBarrelSkip{1.401 * Acts::UnitConstants::degree};
    /// @brief Transition bin: BI-BE and BI-EE
    double transBIToExt{2.3 * Acts::UnitConstants::degree};
    /// @brief Transition bin: BI-EM
    double transBIToEM{1.828 * Acts::UnitConstants::degree};
    /// @brief Transition bin: BI-EI
    double transBIToEI{1.5 * Acts::UnitConstants::degree};
    /// @brief Transition bin: BM-BE, BM-EE, BM-EM and BM-EI
    double transBMToEndcap{1.695 * Acts::UnitConstants::degree};
    /// @brief Transition bin: BO-BE and BO-EE
    double transBOToExt{1.042 * Acts::UnitConstants::degree};
    /// @brief Transition bin: BO-EM
    double transBOToEM{0.772 * Acts::UnitConstants::degree};
    /// @brief Transition bin: BO-EI
    double transBOToEI{1.207 * Acts::UnitConstants::degree};
    /// @brief Transition bin: BE-EE
    double transBEToEE{1.75 * Acts::UnitConstants::degree};
    /// @brief Transition bin: BE-EI and EE-EI
    double transExtToEI{2.3 * Acts::UnitConstants::degree};
    /// @brief Transition bin: BE-EM and EE-EM
    double transExtToEM{0.955 * Acts::UnitConstants::degree};
    /// @brief Transition bin: BE-EO and EE-EO
    double transExtToEO{2.6 * Acts::UnitConstants::degree};
    /// @brief Transition bin: EI-EM and EM-EO
    double transEndcapAdj{1.45 * Acts::UnitConstants::degree};
    /// @brief Transition bin: EI-EO
    double transEndcapSkip{2.6 * Acts::UnitConstants::degree};

    /// @brief Endcap bin: two hits in EI, or two in EO
    double endcapSameInnerOuter{0.9 * Acts::UnitConstants::degree};
    /// @brief Endcap bin: two hits in EM
    double endcapSameMiddle{0.67 * Acts::UnitConstants::degree};
    /// @brief Endcap bin: EI-EM and EM-EO
    double endcapAdjacent{0.607 * Acts::UnitConstants::degree};
    /// @brief Endcap bin: EI-EO
    double endcapSkip{2.6 * Acts::UnitConstants::degree};
  };
  /// @brief Bins of the seed's |eta|, each with its own windows
  enum class EtaBin : std::int8_t { Barrel, Transition, Endcap, nBins };
  /// @brief Half-window per hit station
  using Row_t = std::array<double, PatternTopology::nGroups>;
  /// @brief Half-windows per seed station and hit station
  using Matrix_t = std::array<Row_t, PatternTopology::nGroups>;
  /// @brief Half-windows per |eta| bin, seed station and hit station
  using Windows_t = std::array<Matrix_t, Acts::toUnderlying(EtaBin::nBins)>;
  /// @brief Constructor
  /// @param config Configuration of the seed selection
  explicit SeedSelector(const Config& config);
  /// @brief Whether a hit may seed a pattern: it lies in one of the seeding
  ///        station layers and is not an MDT hit, unless seedFromMdt is set
  /// @param hit Hit of interest
  /// @return True if the hit may seed a pattern
  bool goodForSeeding(const HitPayload& hit) const;
  /// @brief Theta half-windows [rad] to search for hits compatible with a seed
  /// @param seed Seed hit
  /// @return The half-window for the hits of each station, from the table of
  ///         the seed's station and |eta| bin
  const Row_t& thetaSearchWindow(const HitPayload& seed) const;
  /// @brief Bin of a seed's |eta|
  /// @param seed Seed hit
  /// @return The barrel, transition or endcap bin
  EtaBin etaBin(const HitPayload& seed) const;
  /// @brief Largest theta half-window of all seeds and stations
  double maxThetaWindow() const;

 private:
  /// @brief Station layers of the seeds: middle and outer, and inner if
  ///        seedFromInner is set
  std::vector<LayerIndex> m_seedingLayers{LayerIndex::Middle,
                                          LayerIndex::Outer};
  /// @brief Half-windows per |eta| bin, seed station and hit station
  Windows_t m_thetaWindows{};
  /// @brief Configuration of the seed selection
  Config m_cfg;
};
static_assert(Acts::Experimental::detail::PatternSeedSelector<
              SeedSelector, HitPayload, PatternTopology>);

/// @brief Provides the hits that measure only phi, which the global pattern
///        finder adds to its patterns after building them in eta
struct OnlyPhiHitsProvider {
  /// @brief Pattern state of the global pattern finder
  using PatternStateAux = GlobalPatternFinder_t::PatternStateAux;
  /// @brief Phi-only hits, per group
  using PhiHitsPerGroup =
      std::array<std::vector<HitPayload>, PatternTopology::nGroups>;
  /// @brief Phi-only hits that may belong to a pattern: those in the space
  ///        point buckets that hold the pattern's hits, per group
  /// @param pattern Pattern of interest
  /// @param gctx Geometry context
  /// @return The phi-only hits, per group
  PhiHitsPerGroup getPhiOnlyHits(const PatternStateAux& pattern,
                                 const Acts::GeometryContext& gctx) const;
  /// @brief Tracking geometry to look up the measurement surfaces of the
  ///         collected phi-only hits
  const Acts::TrackingGeometry* trackingGeometry{nullptr};
};
static_assert(
    Acts::Experimental::detail::OnlyPhiHitsProvider<
        OnlyPhiHitsProvider, HitPayload, PatternTopology, ExpandedSector>);

/// @brief Transformation from the sector frame, in which the space points are
///        expressed, into the global frame. The space point carries the
///        transform of its own measurement surface into the sector frame, hence
///        the global frame is reached as surfaceToGlobal * surfaceToSector^-1.
/// @param gctx Geometry context
/// @param surface Measurement surface of the space point
/// @param sp Space point of interest
Acts::Transform3 localToGlobalTransform(const Acts::GeometryContext& gctx,
                                        const Acts::Surface& surface,
                                        const MuonSpacePoint& sp);

/// @brief Measurement surface of a space point, looked up by its geometry
///        identifier. The surface defines the local axes of the hit
/// @param trackingGeometry Geometry holding the muon measurement surfaces
/// @param sp Space point of interest
const Acts::Surface& measurementSurface(
    const Acts::TrackingGeometry& trackingGeometry, const MuonSpacePoint& sp);

/// @brief Transformation from the sector frame of a bucket into the global
///        frame. The space points of a bucket are all expressed in the frame
///        of the sector of the first one, hence the transform is determined
///        once per bucket and shared by its hits.
/// @param gctx Geometry context
/// @param trackingGeometry Geometry holding the muon measurement surfaces
/// @param bucket Bucket of interest, which must not be empty
Acts::Transform3 bucketToGlobalTransform(
    const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometry& trackingGeometry,
    const MuonSpacePointBucket& bucket);

/// @brief Result of the global pattern recognition: the space points forming
///        one pattern, grouped by station, together with the pattern's
///        direction, sector and quality estimates
struct MuonGlobalPattern {
  /// @brief Space points of the pattern per station index
  std::array<std::vector<const MuonSpacePoint*>, PatternTopology::nGroups>
      hitsPerStation{};
  /// @brief Average global theta of the pattern
  double theta{0.};
  /// @brief Average global phi of the pattern
  double phi{0.};
  /// @brief Expanded sector index of the pattern
  ExpandedSector::Index_t sector{0};
  /// @brief Number of precision layers
  unsigned nPrecisionLayers{0};
  /// @brief Number of trigger layers
  unsigned nTriggerLayers{0};
  /// @brief Number of phi layers
  unsigned nPhiLayers{0};
  /// @brief Mean over the eta hits of their squared residual divided by the
  ///        acceptance window of the pattern search
  double meanNormResidual2{0.};
};
/// @brief Container of MuonGlobalPattern
using MuonGlobalPatternContainer = std::vector<MuonGlobalPattern>;

/// @brief Algorithm performing global pattern recognition in the muon
///        spectrometer.
///
/// It builds global patterns of precision and non-precision hits from the
/// muon space points of an event. It first builds patterns in eta and then
/// adds compatible phi-only hits to the patterns. The resulting patterns are
/// written into the event store.
class GlobalPatternFinderAlgorithm final : public IAlgorithm {
 public:
  /// @brief Configuration of the algorithm
  struct Config {
    /// @brief Key of the input space point container
    std::string inSpacePoints{};
    /// @brief Key of the output pattern container
    std::string outPatterns{};
    /// @brief Tracking geometry holding the muon measurement surfaces
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry{};

    /// @brief Toggle the utilization of MDT hits to build patterns
    bool useMdtHits{true};
    /// @brief Number of standard deviations to consider for residual acceptance
    double nResidualSigma{3.};
    /// @brief Residual uncertainty [mm] above which a hit has low confidence
    double lowConfidenceResSigma{50. * Acts::UnitConstants::mm};
    /// @brief Maximum angle [rad] of a hit to the pattern's bending plane for
    ///        the hit to be phi compatible with the pattern
    double phiTolerance{3.5 * Acts::UnitConstants::degree};
    /// @brief Minimum number of trigger layers in the bending direction to
    ///        accept a pattern
    unsigned minTriggerLayers{2};
    /// @brief Minimum number of precision layers in the bending direction to
    ///        accept a pattern
    unsigned minPrecisionLayers{8};
    /// @brief Minimum number of phi layers to accept a pattern
    unsigned minPhiLayers{1};
    /// @brief Minimum number of layers for a station to count as a good
    ///        station
    unsigned minStationLayers{4};
    /// @brief Upper cut on the pattern's mean squared normalized residual. A
    ///        large value, e.g. 10, disables the cut
    double meanNormRes2Cut{3.5};
    /// @brief Maximum number of attempts to build a pattern from hits already
    ///        used in existing patterns
    unsigned maxSeedAttempts{3};
    /// @brief Maximum number of measurement layers of a station a pattern
    ///        branch may miss before it is discarded
    unsigned maxMissLayersInStation{3};
    /// @brief Minimum distance [mm] between two hits to define the pattern line
    ///        through them; closer hits use the beam spot instead
    double minHitDistance4Line{200. * Acts::UnitConstants::mm};
    /// @brief Beam spot radius
    double beamSpotRadius{20. * Acts::UnitConstants::cm};
    /// @brief Beam spot length
    double beamSpotLength{1. * Acts::UnitConstants::m};
  };

  /// @brief Constructor
  /// @param cfg Configuration of the algorithm
  /// @param logger Logger of the algorithm; a default one if null
  /// @throws std::invalid_argument if a collection key or the tracking
  ///         geometry is missing
  explicit GlobalPatternFinderAlgorithm(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// @brief Finds the global patterns of an event
  /// @param ctx Context of the event
  /// @return SUCCESS, with the patterns written to the event store
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// @brief Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// @brief Pattern as returned by the global pattern finder
  using PatternResult = GlobalPatternFinder_t::OutputPattern;

  /// @brief Builds the hit payloads of an event and the search tree over them
  /// @param gctx Geometry context
  /// @param spacepoints Space point buckets of the event
  /// @return The hit payloads and the tree that points to them
  SearchTreeData constructTree(
      const Acts::GeometryContext& gctx,
      const MuonSpacePointContainer& spacepoints) const;
  /// @brief Convert a pattern of the finder into a MuonGlobalPattern
  /// @param candidate Pattern to be converted
  /// @return The converted pattern
  MuonGlobalPattern convertToPattern(const PatternResult& candidate) const;
  /// @brief Convert the patterns of the finder into MuonGlobalPatterns
  /// @param candidates Patterns to be converted
  /// @return The converted patterns
  MuonGlobalPatternContainer convertToPattern(
      const std::vector<PatternResult>& candidates) const;

  /// @brief Configuration of the algorithm
  Config m_cfg;

  /// @brief Input space point container
  ReadDataHandle<MuonSpacePointContainer> m_inSpacePoints{this,
                                                          "InSpacePoints"};
  /// @brief Output pattern container
  WriteDataHandle<MuonGlobalPatternContainer> m_outPatterns{this,
                                                            "OutPatterns"};

  /// @brief Seed selector
  std::unique_ptr<SeedSelector> m_seedSelector{};
  /// @brief Provider of the phi-only hits
  OnlyPhiHitsProvider m_onlyPhiProvider{};
  /// @brief Global pattern finder
  std::unique_ptr<GlobalPatternFinder_t> m_globPatFinder{};
};

}  // namespace ActsExamples

ACTS_OSTREAM_FORMATTER(ActsExamples::ExpandedSector::SectorProjector);

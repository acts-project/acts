// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/GlobalPatternFinderAlgorithm.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <format>
#include <iterator>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

using namespace Acts::UnitLiterals;

namespace {

/// @brief Sector geometry of the muon spectrometer used by ExpandedSector:
///        sizes, central phi and overlaps of the 16 sectors
namespace sectorMap {

constexpr double s_oneEightsOfPi{std::numbers::pi / 8.};         //  pi/8
constexpr double s_inverseOneEightsOfPi{8. / std::numbers::pi};  // 8/pi
constexpr std::array<double, 2> s_sectorSize{
    0.4 * s_oneEightsOfPi,
    0.6 * s_oneEightsOfPi};  // side of a sector in radiants
constexpr double s_sectorOverlap{
    0.1 *
    s_oneEightsOfPi};  // size of the overlap between small and large sectors

/// sector size (exclusive) in radians
double sectorSize(int sector) {
  const int idx = sector % 2;
  return s_sectorSize[idx];
}

/// sector width (with overlap) in radians
double sectorWidth(int sector) {
  return sectorSize(sector) + s_sectorOverlap;
}

/// returns the centeral phi position of a sector in radians
double sectorPhi(int sector) {
  if (sector < 10)
    return std::numbers::pi * (sector - 1) / 8.;
  return -std::numbers::pi * (2 - (sector - 1) / 8.);
}

/// transforms a phi position from and to the sector coordinate system in
/// radians
double transformPhiToSector(double phi, int sector, bool toSector = true) {
  double sign = toSector ? -1 : 1;
  double dphi = phi + sign * sectorPhi(sector);
  if (dphi > std::numbers::pi)
    dphi -= 2 * std::numbers::pi;
  if (dphi < -std::numbers::pi)
    dphi += 2 * std::numbers::pi;
  return dphi;
}

/// checks whether the phi position is consistent with sector
bool insideSector(int sector, double phi) {
  double phiInSec = transformPhiToSector(phi, sector);
  if (phiInSec < -sectorWidth(sector))
    return false;
  if (phiInSec > sectorWidth(sector))
    return false;
  return true;
}

/// returns the sector corresponding to the phi position
int getSector(double phi) {
  // remap phi onto sector structure
  double val = (phi + sectorSize(1)) *
               s_inverseOneEightsOfPi;  // convert to value between -8 and 8,
                                        // shift by width first sector
  if (val < 0)
    val += 16;
  int sliceIndex = static_cast<int>(val / 2);  // large/small wedge
  double valueInSlice = val - 2 * sliceIndex;
  int sector = 2 * sliceIndex + 1;
  if (valueInSlice > 1.2)
    ++sector;
  return sector;
}

/// returns the main sector plus neighboring if the phi position is in an
/// overlap region
void getSectors(double phi, std::vector<int>& sectors) {
  int sector = getSector(phi);
  int sectorNext = sector != 16 ? sector + 1 : 1;
  int sectorBefore = sector != 1 ? sector - 1 : 16;
  if (insideSector(sectorBefore, phi))
    sectors.push_back(sectorBefore);
  if (insideSector(sector, phi))
    sectors.push_back(sector);
  if (insideSector(sectorNext, phi))
    sectors.push_back(sectorNext);
}

/// returns the phi position of the overlap between the two sectors (which have
/// to be neighboring) in radians
double sectorOverlapPhi(int sector1, int sector2) {
  if (sector1 == sector2)
    return sectorPhi(sector1);

  int s1 = std::min(sector1, sector2);
  int s2 = std::max(sector1, sector2);
  if (s2 == 16 && s1 == 1) {
    s1 = 16;
    s2 = 1;
  } else if (std::abs(s1 - s2) > 1) {
    return 0;
  }

  double phi1 = sectorPhi(s1);
  double phio1 = phi1 + sectorSize(s1);
  if (phio1 > std::numbers::pi)
    phio1 -= 2 * std::numbers::pi;

  return phio1;
}

}  // namespace sectorMap

/// @brief Number of sectors in the muon spectrometer (MuonStationIndex)
constexpr unsigned numberOfSectors() {
  return 16;
}

using SectorProjector = ActsExamples::ExpandedSector::SectorProjector;
constexpr std::int8_t nExpanded = 2 * numberOfSectors();

inline std::tuple<unsigned, SectorProjector> msSectorAndProj(
    const std::int8_t expandSector) {
  assert(expandSector >= 0 && expandSector <= nExpanded);
  if (expandSector == 0) {
    return std::make_tuple(numberOfSectors(), SectorProjector::center);
  } else if (expandSector == 1) {
    return std::make_tuple(numberOfSectors(), SectorProjector::rightOverlap);
  }
  const int regSector = expandSector / 2;
  const int backConv = expandSector - 2 * regSector;
  switch (backConv) {
    case 1:
      return std::make_tuple(regSector, SectorProjector::rightOverlap);
    case 0:
      return std::make_tuple(regSector, SectorProjector::center);
    case -1:
      return std::make_tuple(regSector, SectorProjector::leftOverlap);
    default:
      throw std::runtime_error(
          std::format("Cannot deduce the sector overlap {}", expandSector));
  }
  // all cases have been addressed
}

}  // namespace

namespace ActsExamples {

Acts::Transform3 localToGlobalTransform(const Acts::GeometryContext& gctx,
                                        const Acts::Surface& surface,
                                        const MuonSpacePoint& sp) {
  return surface.localToGlobalTransform(gctx) *
         sp.toSectorTransform().inverse();
}

const Acts::Surface& measurementSurface(
    const Acts::TrackingGeometry& trackingGeometry, const MuonSpacePoint& sp) {
  const Acts::Surface* surface{trackingGeometry.findSurface(sp.geometryId())};
  if (surface == nullptr) {
    throw std::runtime_error(
        std::format("measurementSurface() No surface for geometry id {}",
                    sp.geometryId().value()));
  }
  return *surface;
}

Acts::Transform3 bucketToGlobalTransform(
    const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometry& trackingGeometry,
    const MuonSpacePointBucket& bucket) {
  assert(!bucket.empty());
  const MuonSpacePoint& first{bucket.front()};
  return localToGlobalTransform(
      gctx, measurementSurface(trackingGeometry, first), first);
}

ExpandedSector::ExpandedSector(const std::int8_t expSector)
    : m_sector{expSector} {}
ExpandedSector::ExpandedSector(const unsigned msSector,
                               const SectorProjector proj) {
  m_sector = (2 * msSector + Acts::toUnderlying(proj)) % nExpanded;
}

ExpandedSector::ExpandedSector(const double phi) {
  std::vector<int> sectors{};
  sectorMap::getSectors(phi, sectors);
  assert(!sectors.empty());
  if (sectors.size() == 1) {
    (*this) = ExpandedSector{static_cast<unsigned>(sectors[0]),
                             SectorProjector::center};
  } else {
    const int dS = (sectors[1] - sectors[0]) % numberOfSectors();
    assert(std::abs(dS) == 1);
    (*this) = ExpandedSector{static_cast<unsigned>(sectors[0]),
                             static_cast<SectorProjector>(dS)};
  }
}

SectorProjector ExpandedSector::projector() const {
  return std::get<1>(msSectorAndProj(sector()));
}
unsigned ExpandedSector::msSector() const {
  return std::get<0>(msSectorAndProj(sector()));
}
std::int8_t ExpandedSector::sector() const {
  return m_sector;
}
unsigned ExpandedSector::adjacentMsSector() const {
  const auto [msSec, proj] = msSectorAndProj(sector());
  if (msSec == 1 && proj == SectorProjector::leftOverlap) {
    return numberOfSectors();
  } else if (msSec == numberOfSectors() &&
             proj == SectorProjector::rightOverlap) {
    return 1;
  }
  return msSec + Acts::toUnderlying(proj);
}
double ExpandedSector::sectorSize() const {
  unsigned sector1{msSector()};
  unsigned sector2{adjacentMsSector()};

  if (sector1 == sector2) {
    return sectorMap::sectorSize(sector1);
  }
  // The overlap size is the same for small and large sectors
  return sectorMap::sectorWidth(sector1) - sectorMap::sectorSize(sector1);
}
bool ExpandedSector::operator<(const ExpandedSector& other) const {
  return sector() < other.sector();
}
bool ExpandedSector::operator==(const ExpandedSector& other) const {
  return sector() == other.sector();
}
bool ExpandedSector::operator!=(const ExpandedSector& other) const {
  return sector() != other.sector();
}
double ExpandedSector::phi() const {
  return sectorMap::sectorOverlapPhi(msSector(), adjacentMsSector());
}
Acts::Vector3 ExpandedSector::normalDir() const {
  return Acts::makeDirectionFromPhiTheta(phi() + 90._degree, 90._degree);
}
bool ExpandedSector::isNeighbour(const ExpandedSector& other) const {
  const int dS = std::abs(other.sector() - sector());
  return std::min(dS, nExpanded - dS) <= 1;
}
bool ExpandedSector::insideSector(const double phi) const {
  if (sectorMap::insideSector(msSector(), phi) &&
      sectorMap::insideSector(adjacentMsSector(), phi)) {
    return true;
  }
  return false;
}

std::string ExpandedSector::toString(const SectorProjector proj) {
  switch (proj) {
    using enum SectorProjector;
    case leftOverlap:
      return "leftOverlap";
    case center:
      return "center";
    case rightOverlap:
      return "rightOverlap";
  }
  return "";
}
std::ostream& ExpandedSector::toString(std::ostream& ostr) const {
  ostr << "Expanded sector: " << static_cast<int>(sector()) << " ->  "
       << projector() << " of sector " << msSector();
  return ostr;
}

StIndex toStationIndex(MuonSpacePoint::MuonId::StationName index) {
  using ChIndex = MuonSpacePoint::MuonId::StationName;
  switch (index) {
    case ChIndex::BIS:
    case ChIndex::BIL:
      return StIndex::BI;
    case ChIndex::BMS:
    case ChIndex::BML:
      return StIndex::BM;
    case ChIndex::BOL:
    case ChIndex::BOS:
      return StIndex::BO;
    case ChIndex::BEE:
      return StIndex::BE;
    case ChIndex::EIL:
    case ChIndex::EIS:
      return StIndex::EI;
    case ChIndex::EML:
    case ChIndex::EMS:
      return StIndex::EM;
    case ChIndex::EOL:
    case ChIndex::EOS:
      return StIndex::EO;
    case ChIndex::EEL:
    case ChIndex::EES:
      return StIndex::EE;
    // No station index for these
    case ChIndex::MaxVal:
    case ChIndex::UnDef:
      break;
  };
  return StIndex::StUnknown;
}
LayerIndex toLayerIndex(StIndex index) {
  switch (index) {
    case StIndex::BI:
    case StIndex::EI:
      return LayerIndex::Inner;
    case StIndex::BM:
    case StIndex::EM:
      return LayerIndex::Middle;
    case StIndex::EO:
    case StIndex::BO:
      return LayerIndex::Outer;
    case StIndex::BE:
      return LayerIndex::BarrelExtended;
    case StIndex::EE:
      return LayerIndex::Extended;
    case StIndex::StUnknown:
    case StIndex::StIndexMax:
      break;
  }
  return LayerIndex::LayerIndexMax;
}
bool isBarrel(const StIndex index) {
  switch (index) {
    case StIndex::BI:
    case StIndex::BM:
    case StIndex::BO:
    case StIndex::BE:
      return true;
    default:
      return false;
  }
}

bool isPrecisionHit(const MuonSpacePoint& hit) {
  using enum MuonSpacePoint::MuonId::TechField;
  // Mdt & micromega measurements are precision hits. An sTgc space point is
  //  one if it stems from a strip, which can only be told apart from a pad if
  //  the space point measures the precision coordinate alone.
  return hit.id().technology() == Mdt || hit.id().technology() == Mm ||
         (hit.id().technology() == sTgc && hit.id().measuresEta() &&
          !hit.id().measuresPhi());
}

SeedSelector::SeedSelector(Config&& config) : m_cfg{config} {
  if (m_cfg.seedFromInner) {
    m_seedinglayers.push_back(LayerIndex::Inner);
  }
}

bool SeedSelector::goodForSeeding(const HitPayload& hit) const {
  const LayerIndex hitLayer{toLayerIndex(hit.station)};
  if (std::ranges::find(m_seedinglayers, hitLayer) == m_seedinglayers.end() ||
      (hit.spacePoint()->isStraw() && !m_cfg.seedFromMdt)) {
    return false;
  }
  return true;
}
double SeedSelector::thetaSearchWindow(const HitPayload& hit) const {
  const LayerIndex hitLayer{toLayerIndex(hit.station)};
  return (hitLayer == LayerIndex::Inner || hitLayer == LayerIndex::Outer)
             ? m_cfg.thetaSearchWindow
             : 0.5 * m_cfg.thetaSearchWindow;
}

using MuonId = MuonSpacePoint::MuonId;

HitPayload::HitPayload(const Acts::GeometryContext& gctx,
                       const MuonSpacePoint* sp,
                       const MuonSpacePointBucket* parentBucket,
                       const Acts::Transform3& localToGlobal,
                       const Acts::Surface* measSurface)
    : position{localToGlobal * sp->localPosition()},
      underlyingSp{sp},
      bucket{parentBucket},
      surface{measSurface} {
  // Only the eta strips of two-dimensional TGC measurements can be tilted
  // against their phi strips; their sensor direction follows the tilt
  if (!sp->id().measuresEta() || !sp->id().measuresPhi() || sp->isStraw() ||
      sp->id().technology() != MuonId::TechField::Tgc) {
    return;
  }
  const auto& surfLinearTrf =
      measSurface->localToGlobalTransform(gctx).linear();
  const Acts::Vector3 phiMeasDir = surfLinearTrf.col(Acts::eY);
  const Acts::Vector3 etaMeasDir =
      localToGlobal.rotation() * sp->toNextSensor();
  const Acts::Vector3 phiSensorDir = surfLinearTrf.col(Acts::eX);
  const double c{etaMeasDir.dot(phiMeasDir)};
  if (std::abs(c) > Acts::s_epsilon) {
    stripAngle = std::atan2(etaMeasDir.dot(phiSensorDir), c);
    nonOrthogonalStrips = true;
  }
}
const Acts::Vector3& HitPayload::globalPosition(
    const Acts::GeometryContext& /*gctx*/) const {
  return position;
}
Acts::Vector3 HitPayload::globalSensorDirection(
    const Acts::GeometryContext& gctx) const {
  const auto& surfLinearTrf = surface->localToGlobalTransform(gctx).linear();

  if (spacePoint()->isStraw()) {
    return surfLinearTrf.col(Acts::eZ);
  } else {
    if (nonOrthogonalStrips) {
      return -std::sin(stripAngle) * surfLinearTrf.col(Acts::eY) +
             std::cos(stripAngle) * surfLinearTrf.col(Acts::eX);
    }
    return surfLinearTrf.col(Acts::eY);
  }
}
Acts::Experimental::detail::StripMeasurementDirections
HitPayload::stripMeasurementDirections(
    const Acts::GeometryContext& gctx) const {
  const auto& surfLinearTrf = surface->localToGlobalTransform(gctx).linear();
  // Phi-only strips measure along the surface y axis
  if (!spacePoint()->id().measuresEta()) {
    return {surfLinearTrf.col(Acts::eY), std::nullopt};
  }
  if (spacePoint()->id().measuresPhi()) {
    const Acts::Vector3 etaMeasDir{
        nonOrthogonalStrips ? Acts::Vector3{globalSensorDirection(gctx).cross(
                                  surfLinearTrf.col(Acts::eZ))}
                            : Acts::Vector3{surfLinearTrf.col(Acts::eX)}};
    return {etaMeasDir, Acts::Vector3{surfLinearTrf.col(Acts::eY)}};
  }
  return {surfLinearTrf.col(Acts::eX), std::nullopt};
}
bool HitPayload::operator==(const HitPayload& other) const {
  return spacePoint() == other.spacePoint();
}
bool HitPayload::overlaps(const HitPayload& other) const {
  // A measurement enters a space point for every non-precision hit it was
  // combined with and for every bucket it was sorted into. Those space points
  // share the surface and the precision coordinate and differ only along the
  // sensor, where the non-precision hit placed them
  return spacePoint()->geometryId() == other.spacePoint()->geometryId() &&
         std::abs(spacePoint()->localPosition().y() -
                  other.spacePoint()->localPosition().y()) < Acts::s_epsilon;
}

bool PatternTopology::layerSorter(const HitPayload& hit1,
                                  const HitPayload& hit2) {
  StIndex st1{hit1.station};
  StIndex st2{hit2.station};
  if (st1 == st2) {
    // Hits in the same spectrometer sector
    if (hit1.spacePoint()->id().sameStation(hit2.spacePoint()->id())) {
      if (hit1.locLayer == hit2.locLayer) {
        return hit1.spacePoint()->localPosition().y() <
               hit2.spacePoint()->localPosition().y();
      }
      return hit1.locLayer < hit2.locLayer;
    }
    // Hits in the same station and different sectors. We can have this case
    // for hits in the overlap region of two adjacent sectors.
    const double delta{isBarrel(st1)
                           ? Acts::VectorHelpers::perp(hit1.position) -
                                 Acts::VectorHelpers::perp(hit2.position)
                           : std::abs(hit1.position.z()) -
                                 std::abs(hit2.position.z())};
    if (std::abs(delta) <= Acts::s_epsilon) {
      return hit1.spacePoint()->localPosition().y() <
             hit2.spacePoint()->localPosition().y();
    }
    return delta < 0.;
  }
  LayerIndex layer1{toLayerIndex(st1)};
  LayerIndex layer2{toLayerIndex(st2)};
  using enum LayerIndex;
  if (layer1 == layer2) {
    // Hit in different stations but same station layer. Expected to happen
    // only for Inner and Middle
    if (layer1 == Middle) {
      // If both hits are in the middle layer, the one in the barrel comes
      // first
      return st1 == StIndex::BM;
    }
    if (layer1 == Inner) {
      // If both hits are in the inner layer, we use the global R, since in
      // large sector BI comes first, while in small sector EI comes first.
      return Acts::VectorHelpers::perp(hit1.position) <
             Acts::VectorHelpers::perp(hit2.position);
    }
    throw std::runtime_error(
        "Unexpected to have two pattern-compatible hits one in BO and the "
        "other in EO.");
  }
  if (layer1 == Inner || layer2 == Inner) {
    // If we have one hit in Inner layer for sure it comes first
    return layer1 == Inner;
  }
  if (layer1 == Outer || layer2 == Outer) {
    // If we have one hit in Outer layer for sure it comes last
    return layer2 == Outer;
  }
  if (layer1 == BarrelExtended || layer2 == BarrelExtended) {
    // If we have one hit in BarrelExtended and the other in the Middle layer,
    // the former comes first
    return layer1 == BarrelExtended;
  }
  // If we have one hit in Extended (EE) layer and the other in the Middle
  // layer, it depends if the latter is endcap or barrel
  if (layer1 == Extended) {
    return st2 == StIndex::EM;
  }
  return st1 == StIndex::BM;
}

bool PatternTopology::sameLayer(const HitPayload& hit1,
                                const HitPayload& hit2) {
  if (hit1.station != hit2.station) {
    return false;
  }
  if (hit1.locLayer == hit2.locLayer &&
      hit1.spacePoint()->id().sameStation(hit2.spacePoint()->id())) {
    return true;
  }
  const double delta{isBarrel(hit1.station)
                         ? Acts::VectorHelpers::perp(hit1.position) -
                               Acts::VectorHelpers::perp(hit2.position)
                         : std::abs(hit1.position.z()) -
                               std::abs(hit2.position.z())};
  return std::abs(delta) <= Acts::s_epsilon;
}
PatternTopology::GroupIdx PatternTopology::groupIndex(const HitPayload& hit) {
  return static_cast<GroupIdx>(Acts::toUnderlying(hit.station));
}

OnlyPhiHitsProvider::PhiHitsPerGroup OnlyPhiHitsProvider::getPhiOnlyHits(
    const PatternStateAux& pattern, const Acts::GeometryContext& gctx) const {
  PhiHitsPerGroup phiOnlyHits{};
  for (const auto& [group, hits] : Acts::enumerate(pattern.hitsPerGroup)) {
    std::vector<const MuonSpacePointBucket*> parentBuckets{};
    for (const auto& hit : hits) {
      const MuonSpacePointBucket* bucket{hit->bucket};
      if (std::ranges::find(parentBuckets, bucket) == parentBuckets.end()) {
        parentBuckets.push_back(bucket);
      }
    }
    for (const MuonSpacePointBucket* bucket : parentBuckets) {
      const Acts::Transform3 localToGlobal{
          bucketToGlobalTransform(gctx, *trackingGeometry, *bucket)};

      for (const MuonSpacePoint& h : *bucket) {
        if (!h.id().measuresEta()) {
          phiOnlyHits[group].emplace_back(
              gctx, &h, bucket, localToGlobal,
              &measurementSurface(*trackingGeometry, h));
        }
      }
    }
  }
  return phiOnlyHits;
}

GlobalPatternFinderAlgorithm::GlobalPatternFinderAlgorithm(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("GlobalPatternFinderAlgorithm", std::move(logger)),
      m_cfg{cfg} {
  if (m_cfg.inSpacePoints.empty()) {
    throw std::invalid_argument(
        "GlobalPatternFinderAlgorithm: Missing space point collection");
  }
  if (m_cfg.outPatterns.empty()) {
    throw std::invalid_argument(
        "GlobalPatternFinderAlgorithm: Missing output pattern collection");
  }
  if (!m_cfg.trackingGeometry) {
    throw std::invalid_argument(
        "GlobalPatternFinderAlgorithm: Missing tracking geometry");
  }
  m_onlyPhiProvider.trackingGeometry = m_cfg.trackingGeometry.get();
  m_inSpacePoints.initialize(m_cfg.inSpacePoints);
  m_outPatterns.initialize(m_cfg.outPatterns);

  SeedSelector::Config selectorCfg{};
  selectorCfg.seedFromMdt = m_cfg.seedFromMdt;
  selectorCfg.thetaSearchWindow = m_cfg.thetaSearchWindow;
  selectorCfg.seedFromInner = m_cfg.seedFromInner;
  m_seedSelector = std::make_unique<SeedSelector>(std::move(selectorCfg));

  GlobalPatternFinder_t::Config patCfg{};
  patCfg.nResidualSigma = m_cfg.nResidualSigma;
  patCfg.lowConfidenceResSigma = m_cfg.lowConfidenceResSigma;
  patCfg.maxThetaSizeOverlap = 2. * m_cfg.thetaSearchWindow;
  patCfg.phiTolerance = m_cfg.phiTolerance;
  patCfg.loosenAcceptanceFactor = 2.5;
  patCfg.minStripEtaLayers = m_cfg.minTriggerLayers;
  patCfg.minPrecisionLayers = m_cfg.minPrecisionLayers;
  patCfg.minPhiLayers = m_cfg.minPhiLayers;
  patCfg.minGroupLayers = m_cfg.minStationLayers;
  patCfg.meanNormRes2Cut = m_cfg.meanNormRes2Cut;
  patCfg.maxSeedAttempts = m_cfg.maxSeedAttempts;
  patCfg.maxMissLayersInGroup = m_cfg.maxMissLayersInStation;
  patCfg.minHitDistance4Line = m_cfg.minHitDistance4Line;

  m_globPatFinder = std::make_unique<GlobalPatternFinder_t>(
      std::move(patCfg), this->logger().clone("GlobalPatternFinder"));

  // Print Configuration
  ACTS_LOG_WITH_LOGGER(
      this->logger(), Acts::Logging::DEBUG,
      "Global Pattern Finder Configuration:\n"
          << " Theta search window [rad]: " << m_cfg.thetaSearchWindow << "\n"
          << " Number of residual standard deviations: " << m_cfg.nResidualSigma
          << "\n"
          << " Low confidence residual sigma [mm]: "
          << m_cfg.lowConfidenceResSigma << "\n"
          << " Max missed layer hits in station: "
          << m_cfg.maxMissLayersInStation << "\n"
          << " Min hit distance for line [mm]: " << m_cfg.minHitDistance4Line
          << "\n"
          << " Phi tolerance [deg]: " << m_cfg.phiTolerance / 1_degree << "\n"
          << " Min trigger layers: " << m_cfg.minTriggerLayers << "\n"
          << " Min precision layers: " << m_cfg.minPrecisionLayers << "\n"
          << " Min phi layers: " << m_cfg.minPhiLayers << "\n"
          << " Min station layers: " << m_cfg.minStationLayers << "\n"
          << " Mean norm residual^2 cut: " << m_cfg.meanNormRes2Cut << "\n"
          << " Seed from inner: " << m_cfg.seedFromInner << "\n"
          << " Use MDT hits: " << m_cfg.useMdtHits << "\n"
          << " Seed from MDT: " << m_cfg.seedFromMdt << "\n"
          << " Max seed attempts: " << m_cfg.maxSeedAttempts << "\n"
          << " Beam spot radius: " << m_cfg.beamSpotRadius << "\n"
          << " Beam spot length: " << m_cfg.beamSpotLength << "\n");
}

ProcessCode GlobalPatternFinderAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const MuonSpacePointContainer& inSpacePoints{m_inSpacePoints(ctx)};
  ACTS_DEBUG("Reading " << inSpacePoints.size()
                        << " SP buckets from collection: "
                        << m_cfg.inSpacePoints);

  const Acts::GeometryContext& gctx{ctx.recoGeoContext};

  SearchTreeData treeData{constructTree(gctx, inSpacePoints)};

  GlobalPatternFinder_t::BeamspotInfo beamSpot{};
  beamSpot.position = Acts::Vector3::Zero();
  beamSpot.radius = m_cfg.beamSpotRadius;
  beamSpot.length = m_cfg.beamSpotLength;

  MuonGlobalPatternContainer patterns{
      convertToPattern(m_globPatFinder->findPatterns(
          gctx, treeData.tree, *m_seedSelector, m_onlyPhiProvider, beamSpot))};

  ACTS_DEBUG("Written " << patterns.size()
                        << " GlobalPatterns into the event store.");
  if (logger().doPrint(Acts::Logging::DEBUG)) {
    for (const MuonGlobalPattern& pat : patterns) {
      std::ostringstream hitsPerStation{};
      for (const auto& [st, hits] : Acts::enumerate(pat.hitsPerStation)) {
        if (!hits.empty()) {
          hitsPerStation << " " << static_cast<int>(st) << ":" << hits.size();
        }
      }
      ACTS_DEBUG("Pattern in "
                 << ExpandedSector{pat.sector} << ", theta: " << pat.theta
                 << ", phi: " << pat.phi
                 << ", precision/trigger/phi layers: " << pat.nPrecisionLayers
                 << "/" << pat.nTriggerLayers << "/" << pat.nPhiLayers
                 << ", meanNormRes2: " << pat.meanNormResidual2
                 << ", hits per station:" << hitsPerStation.str());
    }
  }
  // Write out the global patterns.
  m_outPatterns(ctx, std::move(patterns));
  return ProcessCode::SUCCESS;
}

SearchTreeData GlobalPatternFinderAlgorithm::constructTree(
    const Acts::GeometryContext& gctx,
    const MuonSpacePointContainer& spacepoints) const {
  std::vector<HitPayload> hitPayloads{};
  using enum ExpandedSector::SectorProjector;
  // First estimate the number of hits
  std::size_t totalHits = 0;
  for (const MuonSpacePointBucket& bucket : spacepoints) {
    totalHits += bucket.size();
  }
  hitPayloads.reserve(totalHits);

  for (const MuonSpacePointBucket& bucket : spacepoints) {
    if (bucket.empty()) {
      continue;
    }
    const Acts::Transform3 localToGlobal{
        bucketToGlobalTransform(gctx, *m_cfg.trackingGeometry, bucket)};

    for (const MuonSpacePoint& hit : bucket) {
      // Ignore only-phi hits and MDT hits if desired
      if (!hit.id().measuresEta() || (!m_cfg.useMdtHits && hit.isStraw())) {
        continue;
      }

      const Acts::Surface& surface{
          measurementSurface(*m_cfg.trackingGeometry, hit)};
      hitPayloads.emplace_back(gctx, &hit, &bucket, localToGlobal, &surface);

      if (logger().doPrint(Acts::Logging::VERBOSE)) {
        const HitPayload& newHit{hitPayloads.back()};
        std::ostringstream oss{};
        oss << __func__ << "() Building hit from " << hit << std::endl
            << "Pos: " << Acts::toString(newHit.position) << ", SensorDir: "
            << Acts::toString(newHit.globalSensorDirection(gctx));
        if (!newHit.spacePoint()->isStraw()) {
          oss << ", orthogonalStrips: " << !newHit.nonOrthogonalStrips;
        }
        ACTS_VERBOSE(oss.str());
      }
    }
  }

  SearchTree_t::vector_t treeData{};
  // An event without eta hits gives an empty tree
  if (hitPayloads.empty()) {
    return SearchTreeData{std::move(hitPayloads),
                          SearchTree_t{std::move(treeData)}};
  }
  treeData.reserve(3 * hitPayloads.size());

  // The sector entries of a hit: its own sector and the overlaps with the
  // neighbours, common to all hits of a bucket
  auto addEntry = [&](const double theta, const ExpandedSector& expSect,
                      const HitPayload& hit) {
    std::array<double, 2> coords{};
    using HitCoords = GlobalPatternFinder_t::HitCoords;
    coords[Acts::toUnderlying(HitCoords::eTheta)] = theta;
    coords[Acts::toUnderlying(HitCoords::eSector)] = expSect.sector();
    ACTS_VERBOSE(__func__ << "() Add hit with coordinates [" << coords[0]
                          << ", " << coords[1] << "] to search tree");
    treeData.emplace_back(std::move(coords), &hit);
  };
  for (const HitPayload& hit : hitPayloads) {
    ACTS_VERBOSE(__func__ << "() Spacepoint: " << *hit.spacePoint());
    const Acts::Vector3& pos{hit.position};
    // The sector of the bucket, from its first hit
    const unsigned msSector{hit.bucket->front().id().sector()};
    if (hit.spacePoint()->id().measuresPhi()) {
      // No projection: theta is the same for all sectors
      const double phi{Acts::VectorHelpers::phi(pos)};
      const double theta{Acts::VectorHelpers::theta(pos)};
      for (const ExpandedSector::SectorProjector proj :
           {leftOverlap, center, rightOverlap}) {
        const ExpandedSector expSect{msSector, proj};
        if (proj != center && !expSect.insideSector(phi)) {
          continue;
        }
        addEntry(theta, expSect, hit);
      }
      continue;
    }
    // Project the hit onto the plane along the sector radial direction, which
    // removes the bias of the hit displacement in phi
    const double r2{pos.head<2>().squaredNorm()};
    for (const ExpandedSector::SectorProjector proj :
         {leftOverlap, center, rightOverlap}) {
      const ExpandedSector expSect{msSector, proj};
      const double projR{std::sqrt(
          std::max(0., r2 - Acts::square(pos.dot(expSect.normalDir()))))};
      addEntry(std::atan2(projR, pos.z()), expSect, hit);
    }
  }
  ACTS_VERBOSE(__func__ << "() Create a new tree with " << treeData.size()
                        << " entries and " << hitPayloads.size() << " hits.");
  return SearchTreeData{std::move(hitPayloads),
                        SearchTree_t{std::move(treeData)}};
}

MuonGlobalPattern GlobalPatternFinderAlgorithm::convertToPattern(
    const PatternResult& candidate) const {
  MuonGlobalPattern pattern{};
  // Add eta hits
  for (PatternTopology::GroupIdx g = 0u; g < PatternTopology::nGroups; ++g) {
    const auto& hits{candidate.hitsPerGroup[g]};
    if (hits.empty())
      continue;

    auto& outHits{pattern.hitsPerStation[g]};
    outHits.reserve(hits.size());

    std::ranges::for_each(hits, [&outHits](const HitPayload* h) {
      outHits.push_back(h->spacePoint());
    });
  }

  // Add phi-only hits
  for (const HitPayload& hit : candidate.phiOnlyHits) {
    pattern.hitsPerStation[Acts::toUnderlying(hit.station)].push_back(
        hit.spacePoint());
  }
  pattern.theta = candidate.patTheta;
  pattern.phi = candidate.patPhi;
  // Set the pattern sector(s) and theta.
  pattern.sector = candidate.sector.sector();
  // Set pattern quality information.
  pattern.nPrecisionLayers = candidate.nPrecisionLayers;
  pattern.nTriggerLayers = candidate.nTriggerLayers;
  pattern.nPhiLayers = candidate.nPhiLayers;
  pattern.meanNormResidual2 = candidate.meanNormResidual2;
  return pattern;
}
MuonGlobalPatternContainer GlobalPatternFinderAlgorithm::convertToPattern(
    const std::vector<PatternResult>& candidates) const {
  MuonGlobalPatternContainer patterns{};
  patterns.reserve(candidates.size());
  std::transform(candidates.begin(), candidates.end(),
                 std::back_inserter(patterns),
                 [this](const PatternResult& candidate) {
                   return convertToPattern(candidate);
                 });
  return patterns;
}

}  // namespace ActsExamples

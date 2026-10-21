// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Seeding/GlobalPatternFinder.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <compare>
#include <format>
#include <future>
#include <memory>
#include <mutex>
#include <numbers>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TTree.h>

#include "StrawHitGeneratorHelper.hpp"

/// Unit test of the Global Pattern Finder (GPF) on a toy muon spectrometer
///  with a barrel and two endcaps. In each event, a few muons of fixed energy
///  are propagated from the origin through an ideal toroidal field. In each
///  crossed chamber, the space points are created by the MeasurementGenerator
///  of the straw line fitter tests, using the tangent to the muon trajectory at
///  the station. The GPF runs on the space points of all the muons and the
///  found patterns are compared to the truth. All the GPF ingredients (hit,
///  sector, topology, seed selector and phi-only hit provider) are mocked here,
///  mirroring the ATLAS implementation in
///  MuonFastRecoHelpers/GlobalPatternFinderDefs.h

using TimePoint_t = std::chrono::system_clock::time_point;
using GenCfg_t = ActsTests::MeasurementGenerator::Config;

constexpr auto logLvl = Acts::Logging::Level::INFO;
constexpr auto gpfLogLvl = Acts::Logging::Level::INFO;
/// Number of events per test
constexpr std::size_t nEvents = 10;
/// Number of tests run in parallel
const unsigned nThreads =
    logLvl != Acts::Logging::Level::INFO
        ? 1u
        : std::max(std::thread::hardware_concurrency(), 1u);
std::mutex writeMutex{};

ACTS_LOCAL_LOGGER(getDefaultLogger("GlobalPatternFinderTest", logLvl));

namespace ActsTests {
using namespace Acts::UnitLiterals;
///////////////////////////////////////////////////////////////////////////////
/// Toy detector description
///////////////////////////////////////////////////////////////////////////////
/// @brief Muon stations. Each station is one group of the pattern topology.
///         There is no transition region: barrel and endcap muons cross only
///         the barrel and only the endcap stations, respectively.
enum class Station : std::uint8_t { BI = 0, BM, BO, EI, EM, EO, nStations };
constexpr std::size_t s_nStations = toUnderlying(Station::nStations);

constexpr bool isBarrel(const Station st) {
  return st == Station::BI || st == Station::BM || st == Station::BO;
}
/// @brief Station layer: 0 = Inner, 1 = Middle, 2 = Outer
constexpr std::uint8_t stationLayer(const Station st) {
  return toUnderlying(st) % 3;
}

/// @brief Azimuthal sectorization: 2pi is divided into 16 equal sectors,
///         without overlaps. Sector 0 is centred at phi = 0.
class TestSector {
 public:
  using Index_t = std::uint8_t;
  static constexpr Index_t s_nSectors{16};
  static constexpr double s_sectorSize{360._degree / s_nSectors};

  /// @brief Constructor from the sector index
  explicit TestSector(const Index_t idx)
      : m_idx{static_cast<Index_t>(idx % s_nSectors)} {}
  /// @brief Constructor of the sector containing a given phi.
  explicit TestSector(const double phi)
      : m_idx{static_cast<Index_t>(
            std::lround(Acts::detail::radian_pos(phi) / s_sectorSize) %
            s_nSectors)} {}

  Index_t sector() const { return m_idx; }
  double sectorSize() const { return s_sectorSize; }
  /// @brief Phi of the sector centre, in [-pi, pi)
  double phi() const {
    return Acts::detail::radian_sym(sector() * sectorSize());
  }
  bool isNeighbour(const TestSector& other) const {
    const int diff{(sector() + s_nSectors - other.sector()) % s_nSectors};
    return diff <= 1 || diff == s_nSectors - 1;
  }
  bool insideSector(const double testPhi) const {
    return std::abs(deltaPhi(testPhi)) <= 0.5 * sectorSize();
  }
  /// @brief Signed phi distance from the sector centre, in [-pi, pi]
  double deltaPhi(const double testPhi) const {
    return Acts::detail::difference_periodic(testPhi, phi(), 360._degree);
  }
  /// @brief Direction from the beam axis to the sector centre
  Vector3 radialDir() const {
    return Acts::makeDirectionFromPhiTheta(phi(), 90._degree);
  }
  /// @brief Direction orthogonal to the sector plane (tube direction)
  Vector3 normalDir() const {
    return Acts::makeDirectionFromPhiTheta(phi() + 90._degree, 90._degree);
  }

 private:
  Index_t m_idx{0};
};
static_assert(SectorType<TestSector>);

///////////////////////////////////////////////////////////////////////////////
/// Mocked GPF ingredients
///////////////////////////////////////////////////////////////////////////////
/// @brief GPF hit. It wraps a space point created in the local chamber frame by
///         the MeasurementGenerator and caches its global representation.
struct TestHit {
  const FitTestSpacePoint* sp{nullptr};
  /// @brief Global position
  Vector3 position{Vector3::Zero()};
  /// @brief Global sensor direction (strip in the bending direction for
  ///         combined strip space points)
  Vector3 sensorDir{Vector3::UnitX()};
  /// @brief Global measured direction of the strip
  Vector3 measDir{Vector3::UnitY()};
  Station station{Station::BI};
  TestSector::Index_t sector{0};
  /// @brief Distance of the measurement layer plane from the beam axis
  ///         (barrel) or from the transverse plane through the origin (endcap)
  double layerDistance{0.};
  /// @brief Index of the generating muon (-1 for background)
  int muonIdx{-1};

  const FitTestSpacePoint& spacePoint() const { return *sp; }
  const Vector3& globalPosition(const GeometryContext& /*gctx*/) const {
    return position;
  }
  const Vector3& globalSensorDirection(const GeometryContext& /*gctx*/) const {
    return sensorDir;
  }
  /// The strips are orthogonal: only the primary measured direction is set
  StripMeasurementDirections stripMeasurementDirections(
      const GeometryContext& /*gctx*/) const {
    return {measDir, std::nullopt};
  }
  bool isPrecision() const { return sp->isStraw(); }
  bool overlaps(const TestHit& other) const { return sp == other.sp; }
};
static_assert(GlobPatFinderHit<TestHit>);

/// @brief Pattern topology: one group per station. The stations are ordered
///         Inner -> Middle -> Outer and, inside a station, the layers by their
///         distance from the interaction point.
struct TestTopology {
  static constexpr std::uint8_t nGroups{static_cast<std::uint8_t>(s_nStations)};

  static std::weak_ordering compareLayers(const TestHit& hit1,
                                          const TestHit& hit2) {
    if (hit1.station != hit2.station) {
      if (stationLayer(hit1.station) != stationLayer(hit2.station)) {
        return stationLayer(hit1.station) <=> stationLayer(hit2.station);
      }
      return toUnderlying(hit1.station) <=> toUnderlying(hit2.station);
    }
    if (std::abs(hit1.layerDistance - hit2.layerDistance) < s_epsilon) {
      return std::weak_ordering::equivalent;
    }
    return hit1.layerDistance < hit2.layerDistance
               ? std::weak_ordering::less
               : std::weak_ordering::greater;
  }
  static std::uint8_t groupIndex(const TestHit& hit) {
    return static_cast<std::uint8_t>(toUnderlying(hit.station));
  }
};
static_assert(PatternTopology<TestTopology, TestHit>);

using GlobalPatternFinder_t =
    GlobalPatternFinder<TestHit, TestSector, TestTopology>;
using SearchTree_t = GlobalPatternFinder_t::SearchTree_t;
using HitCoords = GlobalPatternFinder_t::HitCoords;

/// @brief Seed selector: seeds from middle & outer strip hits. A single,
///         large theta window is used for all the station pairs, since the
///         test probes the pattern building and not the CPU performance.
class TestSeedSelector {
 public:
  struct Config {
    bool seedFromMdt{false};
    /// @brief Theta search window, used for all the station pairs
    double thetaWindow{3._degree};
  };
  using Window_t = std::array<double, TestTopology::nGroups>;

  explicit TestSeedSelector(const Config& cfg)
      : m_cfg{cfg},
        m_window{filledArray<double, TestTopology::nGroups>(cfg.thetaWindow)} {}
  bool goodForSeeding(const TestHit& hit) const {
    return stationLayer(hit.station) != 0 &&
           (!hit.sp->isStraw() || m_cfg.seedFromMdt);
  }
  const Window_t& thetaSearchWindow(const TestHit& /*hit*/) const {
    return m_window;
  }
  double maxThetaWindow() const { return m_cfg.thetaWindow; }

 private:
  Config m_cfg{};
  Window_t m_window{};
};
static_assert(PatternSeedSelector<TestSeedSelector, TestHit, TestTopology>);

/// @brief Provides the phi-only hits from the chambers (station & sector)
///         contributing eta hits to the pattern
struct TestPhiProvider {
  using PatternState_t = GlobalPatternFinder_t::PatternStateAux;
  using PhiHitsPerGroup =
      std::array<std::vector<TestHit>, TestTopology::nGroups>;

  const std::vector<TestHit>* phiOnlyHits{nullptr};

  PhiHitsPerGroup getPhiOnlyHits(const PatternState_t& pattern,
                                 const GeometryContext& /*gctx*/) const {
    PhiHitsPerGroup result{};
    for (const auto [group, hits] : enumerate(pattern.hitsPerGroup)) {
      std::vector<TestSector::Index_t> sectors{};
      for (const auto& hit : hits) {
        if (std::ranges::find(sectors, hit->sector) == sectors.end()) {
          sectors.push_back(hit->sector);
        }
      }
      for (const TestHit& phiHit : *phiOnlyHits) {
        if (TestTopology::groupIndex(phiHit) == group &&
            std::ranges::find(sectors, phiHit.sector) != sectors.end()) {
          result[group].push_back(phiHit);
        }
      }
    }
    return result;
  }
};
static_assert(
    OnlyPhiHitsProvider<TestPhiProvider, TestHit, TestTopology, TestSector>);

///////////////////////////////////////////////////////////////////////////////
/// Detector layout & muon generator
///////////////////////////////////////////////////////////////////////////////
/// @brief Position of the first tube layer in the local frame of the
///         MeasurementGenerator (chamberDistance in its spawn method)
constexpr double s_genFirstTubeZ{5._m};
/// @brief Local z of the station centre in the MeasurementGenerator frame,
///         half way between the first and the last tube layer
constexpr double s_genStationCentreZ{5.8_m};

/// @brief Configuration of the hit generation in each chamber: the
///         MeasurementGenerator configuration, extended with the options of
///         this test. It is passed to the MeasurementGenerator as its base
///         class.
struct HitGeneratorConfig : GenCfg_t {
  /// @brief MDT + combined trigger strips (RPC/TGC like): 2 gas gaps below
  ///         and above the tubes
  HitGeneratorConfig() {
    createStraws = true;
    twinStraw = false;
    createStrips = true;
    combineSpacePoints = true;
    smearStrips = false;
    discretizeStrips = true;
    stripPitchLoc0 = 3._cm;
    stripPitchLoc1 = 3._cm;
    stripDirLoc0 = std::vector<Vector3>(2, Vector3::UnitY());
    stripDirLoc1 = std::vector<Vector3>(2, Vector3::UnitX());
  }
  /// @brief Number of phi-only strip layers on each side of the tubes (phi
  ///         strips without an eta partner), outside the combined strip
  ///         layers
  std::size_t nPhiOnlyLayers{1};
  /// @brief Add background: in a few stations crossed by each muon, an
  ///         independent straight line crossing the muon trajectory
  bool addBackground{false};
  /// @brief Number of stations per muon with a background line
  std::size_t nBkgStationsPerMuon{2};
  /// @brief Minimum and maximum angle between a background line and the
  ///         muon in the bending plane
  double bkgMinAngle{5._degree};
  double bkgMaxAngle{40._degree};
  /// @brief Maximum distance of a background line from the muon in the
  ///         crossed station
  double bkgMaxOffset{30._cm};
  /// @brief Maximum inclination of a background line in the non-bending
  ///         plane
  double bkgMaxTilt{5._degree};
};

struct DetectorConfig {
  /// @brief Distance of the first tube layer of each station from the beam
  ///         axis (barrel) or from the transverse plane through the origin
  ///         (endcap). Each station is about 2.3 m thick.
  std::array<double, s_nStations> firstTubeDist{4.65_m, 7.0_m,  9.3_m,
                                                7.0_m,  13.0_m, 20.0_m};
  /// @brief Station extension along |z| (barrel) or along rho (endcap)
  std::array<double, s_nStations> extMin{0._m, 0._m, 0._m, 1._m, 1._m, 1._m};
  std::array<double, s_nStations> extMax{6.5_m, 9.5_m, 12.5_m,
                                         5._m,  12._m, 12.5_m};
  /// @brief Barrel toroid: rhoMin < rho < rhoMax and |z| < endcap zMin
  double barrelRhoMin{5.1_m};
  double barrelRhoMax{10.5_m};
  double barrelField{0.5_T};
  /// @brief Endcap toroids: zMin < |z| < zMax and rho < barrel rhoMax
  double endcapZMin{9.2_m};
  double endcapZMax{12.4_m};
  double endcapField{0.5_T};
  /// @brief Step size of the trajectory stepping and limits at which it
  ///         stops
  double trajStepSize{5._mm};
  double trajMaxRho{14._m};
  double trajMaxAbsZ{24._m};
  double trajMaxPath{40._m};
  /// @brief Configuration of the hit generation, used for each chamber
  HitGeneratorConfig genCfg{};

  /// @brief Field strength (along phi) at a given position
  double field(const double rho, const double z) const {
    const double absZ{std::abs(z)};
    if (rho > barrelRhoMax) {
      return 0.;
    }
    if (absZ >= endcapZMin) {
      return absZ < endcapZMax ? endcapField : 0.;
    }
    return rho > barrelRhoMin ? barrelField : 0.;
  }
};

enum class Region : std::uint8_t { Barrel, Endcap };

/// @brief Truth information of the generated muon
struct MuonTruth {
  double energy{0.};
  double p{0.};
  double eta{0.};
  double phi{0.};
  int charge{1};
  double theta() const { return 2. * std::atan(std::exp(-eta)); }
};

struct GeneratorConfig {
  /// @brief Muon energy
  double energy{20._GeV};
  /// @brief |eta| range of the muons in the barrel and in the endcap
  std::array<double, 2> barrelEta{0., 0.9};
  std::array<double, 2> endcapEta{1.5, 2.3};
};

/// @brief Random sign, -1 or +1 with equal probability
double randomSign(RandomEngine& engine) {
  return uniform{0., 1.}(engine) < 0.5 ? -1. : 1.;
}

MuonTruth generateMuon(RandomEngine& engine, const GeneratorConfig& cfg,
                       const Region region) {
  constexpr double muonMass{105.6583755_MeV};
  MuonTruth muon{};
  muon.energy = cfg.energy;
  muon.p = std::sqrt(square(cfg.energy) - square(muonMass));
  const auto& etaRange{region == Region::Barrel ? cfg.barrelEta
                                                : cfg.endcapEta};
  const double absEta{uniform{etaRange[0], etaRange[1]}(engine)};
  muon.eta = randomSign(engine) * absEta;
  muon.phi = uniform{-std::numbers::pi, std::numbers::pi}(engine);
  muon.charge = static_cast<int>(randomSign(engine));
  return muon;
}

/// @brief Trajectory of a muon from the origin in the half plane at constant
///         phi. The toroidal field is along phi, so the muon bends in the
///         (rho, z) plane and phi stays constant. The direction is parametrized
///         by the angle alpha w.r.t. the transverse plane. Along the path s:
///         d(alpha)/ds = q * B / p, d(rho)/ds = cos(alpha), dz/ds = sin(alpha).
///         The trajectory is stepped with exact arcs of constant curvature.
class ToroidTrajectory {
 public:
  struct Point {
    double rho{0.};
    double z{0.};
    double alpha{0.};
  };

  ToroidTrajectory(const MuonTruth& muon, const DetectorConfig& cfg) {
    const double stepSize{cfg.trajStepSize};
    Point point{0., 0., std::atan(std::sinh(muon.eta))};
    m_points.push_back(point);
    for (double s = 0.; s < cfg.trajMaxPath; s += stepSize) {
      const double kappa{muon.charge * cfg.field(point.rho, point.z) / muon.p};
      const double alpha1{point.alpha + kappa * stepSize};
      if (std::abs(kappa) < s_epsilon) {
        point.rho += stepSize * std::cos(point.alpha);
        point.z += stepSize * std::sin(point.alpha);
      } else {
        point.rho += (std::sin(alpha1) - std::sin(point.alpha)) / kappa;
        point.z -= (std::cos(alpha1) - std::cos(point.alpha)) / kappa;
      }
      point.alpha = alpha1;
      /// Stop if the muon leaves the spectrometer or curls back
      if (point.rho < 0. || point.rho > cfg.trajMaxRho ||
          std::abs(point.z) > cfg.trajMaxAbsZ || std::cos(point.alpha) < 0.) {
        break;
      }
      m_points.push_back(point);
    }
  }
  /// @brief Return the first point where the function f of the trajectory
  ///         points changes sign from negative to positive. The point is
  ///         linearly interpolated between the two closest steps.
  template <typename Func_t>
  std::optional<Point> crossing(const Func_t& f) const {
    for (std::size_t i = 1; i < m_points.size(); ++i) {
      const double f0{f(m_points[i - 1])};
      const double f1{f(m_points[i])};
      if (f0 < 0. && f1 >= 0.) {
        const double w{f0 / (f0 - f1)};
        const Point& p0{m_points[i - 1]};
        const Point& p1{m_points[i]};
        return Point{p0.rho + w * (p1.rho - p0.rho), p0.z + w * (p1.z - p0.z),
                     p0.alpha + w * (p1.alpha - p0.alpha)};
      }
    }
    return std::nullopt;
  }

 private:
  std::vector<Point> m_points{};
};

/// @brief Local frame of a chamber of the MeasurementGenerator: x along the
///         tubes, y along the bending direction and z orthogonal to the tube
///         layers, pointing outwards. The station plane at the local
///         z = s_genFirstTubeZ is placed at the first tube layer distance.
///         Barrel: y = global z, z = sector radial direction.
///         Endcap: y = sector radial direction, z = +- global z.
/// @brief Offset between the local z of the MeasurementGenerator and the
///         distance from the beam axis (barrel) or from the transverse plane
///         through the origin (endcap)
double chamberOffset(const Station st, const DetectorConfig& cfg) {
  return cfg.firstTubeDist[toUnderlying(st)] - s_genFirstTubeZ;
}

Transform3 chamberFrame(const Station st, const TestSector& sector,
                        const double zSign, const DetectorConfig& cfg) {
  const double offset{chamberOffset(st, cfg)};
  Transform3 localToGlobal{Transform3::Identity()};
  if (isBarrel(st)) {
    localToGlobal.linear().col(ePos0) = sector.normalDir();
    localToGlobal.linear().col(ePos1) = Vector3::UnitZ();
    localToGlobal.linear().col(ePos2) = sector.radialDir();
    localToGlobal.translation() = offset * sector.radialDir();
  } else {
    const Vector3 normal{zSign * Vector3::UnitZ()};
    localToGlobal.linear().col(ePos0) = sector.radialDir().cross(normal);
    localToGlobal.linear().col(ePos1) = sector.radialDir();
    localToGlobal.linear().col(ePos2) = normal;
    localToGlobal.translation() = offset * normal;
  }
  return localToGlobal;
}

/// @brief Create the space points of a straight line crossing a station. The
///         MeasurementGenerator is run on the line expressed in the local
///         chamber frame.
///  @param globPos: Point on the line, in global coordinates
///  @param globDir: Direction of the line, in global coordinates
///  @param muonIdx: Index of the generating muon (-1 for background)
///  @param outSpacePoints: Container owning the generated space points
///  @param outHits: Hits wrapping the generated space points
void spawnInStation(const Vector3& globPos, const Vector3& globDir,
                    const Station station, const TestSector& sector,
                    const double zSign, const int muonIdx,
                    const DetectorConfig& cfg, RandomEngine& engine,
                    Container_t& outSpacePoints,
                    std::vector<TestHit>& outHits) {
  const auto st{toUnderlying(station)};
  const Transform3 localToGlobal{chamberFrame(station, sector, zSign, cfg)};
  const Vector3 locPos{localToGlobal.inverse() * globPos};
  const Vector3 locDir{localToGlobal.linear().transpose() * globDir};

  /// Line parameters are defined at the local z = 0 plane
  Line_t::ParamVector linePars{};
  using ParIndex = Line_t::ParIndex;
  linePars[toUnderlying(ParIndex::x0)] =
      locPos.x() - locPos.z() * locDir.x() / locDir.z();
  linePars[toUnderlying(ParIndex::y0)] =
      locPos.y() - locPos.z() * locDir.y() / locDir.z();
  linePars[toUnderlying(ParIndex::theta)] = VectorHelpers::theta(locDir);
  linePars[toUnderlying(ParIndex::phi)] = VectorHelpers::phi(locDir);
  const Line_t line{linePars};

  /// The MeasurementGenerator only takes its own configuration: the phi-only
  /// strips are the loc0 strip layers beyond the loc1 ones
  GenCfg_t genCfg{cfg.genCfg};
  genCfg.stripDirLoc0.resize(
      genCfg.stripDirLoc1.size() + cfg.genCfg.nPhiOnlyLayers, Vector3::UnitY());
  genCfg.localToGlobal = localToGlobal;
  for (auto& sp :
       MeasurementGenerator::spawn(line, 0., engine, genCfg, logger())) {
    TestHit hit{};
    hit.sp = sp.get();
    /// Phi-only strips are created at the local y = 0: the strip covers the
    /// whole chamber along the bending direction, so place it at the chamber
    /// centre and store its half length. For barrel chambers, the chamber is
    /// in the z hemisphere crossed by the track.
    const bool isPhiOnly{!sp->isStraw() && !sp->measuresLoc1()};
    const Vector3 locTrackPos{locPos + (sp->localPosition().z() - locPos.z()) /
                                           locDir.z() * locDir};
    if (isPhiOnly) {
      const double ySign{isBarrel(station) && locTrackPos.y() < 0. ? -1. : 1.};
      sp->updatePosition(
          Vector3{sp->localPosition().x(),
                  ySign * 0.5 * (cfg.extMin[st] + cfg.extMax[st]),
                  sp->localPosition().z()});
      sp->setStripHalfLength(0.5 * (cfg.extMax[st] - cfg.extMin[st]));
    }
    hit.position = localToGlobal * sp->localPosition();
    /// The generated chambers are unbounded: apply the station extension. For
    /// phi-only strips, the track position on the strip plane is checked.
    const Vector3 extPos{isPhiOnly ? Vector3{localToGlobal * locTrackPos}
                                   : hit.position};
    const double ext{isBarrel(station) ? std::abs(extPos.z())
                                       : VectorHelpers::perp(extPos)};
    if (ext < cfg.extMin[st] || ext > cfg.extMax[st]) {
      continue;
    }
    /// Combined strip space points are created with (loc0, loc1) strip
    /// directions: the strip in the bending direction is the second one
    const Vector3 locSensorDir{!sp->isStraw() && sp->measuresLoc0() &&
                                       sp->measuresLoc1()
                                   ? sp->toNextSensor()
                                   : sp->sensorDirection()};
    hit.sensorDir = localToGlobal.linear() * locSensorDir;
    hit.measDir = localToGlobal.linear() * locSensorDir.cross(Vector3::UnitZ());
    hit.station = station;
    hit.sector = sector.sector();
    hit.layerDistance = chamberOffset(station, cfg) + sp->localPosition().z();
    hit.muonIdx = muonIdx;
    outHits.push_back(hit);
    outSpacePoints.push_back(std::move(sp));
  }
}

/// @brief Distance of the centre plane of a station from the beam axis
///         (barrel) or from the transverse plane through the origin (endcap)
double stationCentreDist(const Station station, const DetectorConfig& cfg) {
  return chamberOffset(station, cfg) + s_genStationCentreZ;
}

/// @brief Point and tangent of the muon trajectory at the centre plane of each
///         station of its sector (nullopt if the station is not crossed)
using StationCrossings =
    std::array<std::optional<std::pair<Vector3, Vector3>>, s_nStations>;
StationCrossings muonCrossings(const MuonTruth& muon,
                               const DetectorConfig& cfg) {
  const ToroidTrajectory trajectory{muon, cfg};
  const Vector3 muonRadialDir{std::cos(muon.phi), std::sin(muon.phi), 0.};
  const TestSector sector{muon.phi};
  const double cosDPhi{std::cos(sector.deltaPhi(muon.phi))};
  StationCrossings crossings{};
  for (std::size_t st = 0; st < s_nStations; ++st) {
    const auto station{static_cast<Station>(st)};
    /// The muon moves in the half plane at its own phi, which crosses the
    /// barrel planes of the sector at rho * cos(dPhi) = centreDist.
    const double centreDist{stationCentreDist(station, cfg)};
    const auto point{trajectory.crossing([&](const ToroidTrajectory::Point& p) {
      return isBarrel(station) ? p.rho * cosDPhi - centreDist
                               : std::abs(p.z) - centreDist;
    })};
    if (point) {
      crossings[st] = std::make_pair(
          point->rho * muonRadialDir + point->z * Vector3::UnitZ(),
          std::cos(point->alpha) * muonRadialDir +
              std::sin(point->alpha) * Vector3::UnitZ());
    }
  }
  return crossings;
}

/// @brief Side of the muon in z: -1 or +1
double zSide(const MuonTruth& muon) {
  return muon.eta < 0. ? -1. : 1.;
}

/// @brief Create the space points of the muon in the chambers of its sector.
///         In each station, the MeasurementGenerator is run on the tangent to
///         the trajectory at the station centre plane.
///  @param crossings: Crossings of the muon with the station centre planes
///  @param outSpacePoints: Container owning the generated space points
///  @return The hits wrapping the generated space points
std::vector<TestHit> spawnHits(const MuonTruth& muon,
                               const StationCrossings& crossings,
                               const int muonIdx, const DetectorConfig& cfg,
                               RandomEngine& engine,
                               Container_t& outSpacePoints) {
  const TestSector sector{muon.phi};
  const double zSign{zSide(muon)};
  std::vector<TestHit> hits{};
  for (const auto [st, crossing] : enumerate(crossings)) {
    if (crossing) {
      spawnInStation(crossing->first, crossing->second,
                     static_cast<Station>(st), sector, zSign, muonIdx, cfg,
                     engine, outSpacePoints, hits);
    }
  }
  return hits;
}

/// @brief Create the background space points crossing the pattern of the
///         muon. A few random stations crossed by the muon are picked and, in
///         each of them, an independent straight line is generated: different
///         stations get uncorrelated lines. Each line passes close to the point
///         where the muon crosses the station, with a random inclination w.r.t.
///         the muon in the bending plane and a small one in the non-bending
///         plane.
///  @param crossings: Crossings of the muon with the station centre planes
///  @param outSpacePoints: Container owning the generated space points
///  @return The hits wrapping the generated space points
std::vector<TestHit> spawnBackground(const MuonTruth& muon,
                                     const StationCrossings& crossings,
                                     const DetectorConfig& cfg,
                                     RandomEngine& engine,
                                     Container_t& outSpacePoints) {
  const HitGeneratorConfig& genCfg{cfg.genCfg};
  const TestSector sector{muon.phi};
  const double zSign{zSide(muon)};
  std::vector<std::size_t> crossedStations{};
  for (const auto [st, crossing] : enumerate(crossings)) {
    if (crossing) {
      crossedStations.push_back(st);
    }
  }
  std::ranges::shuffle(crossedStations, engine);
  crossedStations.resize(
      std::min(crossedStations.size(), genCfg.nBkgStationsPerMuon));

  std::vector<TestHit> hits{};
  for (const std::size_t st : crossedStations) {
    const auto& [muonPos, muonDir] = *crossings[st];
    /// Directions in the bending plane orthogonal to the muon and tangent to
    /// the sector
    const Vector3 tangentDir{Vector3::UnitZ()
                                 .cross(Vector3{muonPos.x(), muonPos.y(), 0.})
                                 .normalized()};
    const Vector3 bendDir{tangentDir.cross(muonDir)};

    const double angle{randomSign(engine) *
                       uniform{genCfg.bkgMinAngle, genCfg.bkgMaxAngle}(engine)};
    const Vector3 linePos{
        muonPos +
        uniform{-genCfg.bkgMaxOffset, genCfg.bkgMaxOffset}(engine)*bendDir};
    const Vector3 lineDir{
        (std::cos(angle) * muonDir + std::sin(angle) * bendDir +
         std::tan(uniform{-genCfg.bkgMaxTilt, genCfg.bkgMaxTilt}(engine)) *
             tangentDir)
            .normalized()};
    spawnInStation(linePos, lineDir, static_cast<Station>(st), sector, zSign,
                   -1, cfg, engine, outSpacePoints, hits);
  }
  return hits;
}

/// @brief Fill the search tree. Since the sectors do not overlap, each hit is
///         inserted only in its own sector.
SearchTree_t buildSearchTree(const std::vector<TestHit>& etaHits) {
  constexpr auto thetaIdx{toUnderlying(HitCoords::eTheta)};
  constexpr auto sectorIdx{toUnderlying(HitCoords::eSector)};
  SearchTree_t::vector_t entries{};
  entries.reserve(etaHits.size());
  for (const TestHit& hit : etaHits) {
    SearchTree_t::coordinate_t coords{};
    coords[thetaIdx] = VectorHelpers::theta(hit.position);
    coords[sectorIdx] = hit.sector;
    entries.emplace_back(coords, &hit);
  }
  return SearchTree_t{std::move(entries)};
}

///////////////////////////////////////////////////////////////////////////////
/// The test
///////////////////////////////////////////////////////////////////////////////
struct TestConfig {
  std::string name{};
  /// @brief Region of each muon of the event
  std::vector<Region> muonRegions{Region::Barrel, Region::Barrel};
  GeneratorConfig genCfg{};
  DetectorConfig detCfg{};
  TestSeedSelector::Config seedCfg{};
  /// @brief Beamspot passed to the GPF
  GlobalPatternFinder_t::BeamspotInfo beamSpot{Vector3::Zero(), 20._cm, 1._m};
  /// @brief A pattern is matched to the muon contributing at least this
  ///         fraction of its hits
  double patternMatchFraction{0.5};
  /// @brief A muon is found if the hit efficiency and the purity of its best
  ///         pattern are at least these values
  double minFoundEfficiency{0.5};
  double minFoundPurity{0.5};
  unsigned seed{0};
};

struct TestResult {
  /// @brief Fraction of muons whose pattern was found
  double foundFraction{0.};
  /// @brief Number of patterns not matched to any muon
  std::size_t nFakes{0};
  /// @brief Wall-clock time spent in the test [ms]
  long int time{0};
  /// @brief Mean time spent in the pattern finding per event [us]
  double gpfTime{0.};
};

GlobalPatternFinder_t::Config makeFinderConfig(
    const TestSeedSelector& selector) {
  /// Default values of the MuonGlobalPatternFindingAlg
  GlobalPatternFinder_t::Config cfg{};
  cfg.nResidualSigma = 3.;
  cfg.lowConfidenceResSigma = 50.;
  cfg.maxThetaSizeOverlap = 2. * selector.maxThetaWindow() + 0.5_degree;
  cfg.phiTolerance = 3.5_degree;
  cfg.loosenAcceptanceFactor = 2.5;
  cfg.minStripEtaLayers = 2;
  cfg.minPrecisionLayers = 8;
  cfg.minPhiLayers = 1;
  cfg.minGroupLayers = 4;
  cfg.meanNormRes2Cut = 3.5;
  cfg.maxSeedAttempts = 3;
  cfg.maxMissLayersInGroup = 3;
  cfg.minHitDistance4Line = 200._mm;
  return cfg;
}

#define DECLARE_BRANCH(dType, bName) \
  dType bName{};                     \
  outTree->Branch(#bName, &bName);
// NOLINTBEGIN
/// @brief Run the pattern finding on nEvents events. The tree is filled once
///         per generated muon.
TestResult runPatternTest(const TestConfig& testCfg, TFile& outFile) {
  // NOLINTEND
  auto outTree = std::make_unique<TTree>(
      std::format("{:}Tree", testCfg.name).c_str(), "MonitorTree");
  outTree->SetDirectory(nullptr);

  DECLARE_BRANCH(unsigned, event);
  DECLARE_BRANCH(unsigned, nPatterns);
  DECLARE_BRANCH(unsigned, nFakePatterns);
  DECLARE_BRANCH(unsigned, nBkgEtaHits);
  DECLARE_BRANCH(double, gpfTime);

  DECLARE_BRANCH(int, muonIdx);
  DECLARE_BRANCH(char, isBarrelMuon);
  DECLARE_BRANCH(double, trueE);
  DECLARE_BRANCH(double, trueEta);
  DECLARE_BRANCH(double, truePhi);
  DECLARE_BRANCH(double, trueTheta);
  DECLARE_BRANCH(int, trueQ);
  DECLARE_BRANCH(int, trueSector);
  DECLARE_BRANCH(unsigned, nTrueEtaHits);
  DECLARE_BRANCH(unsigned, nTruePhiOnlyHits);
  DECLARE_BRANCH(unsigned, nTrueStations);

  DECLARE_BRANCH(char, found);
  DECLARE_BRANCH(unsigned, nMatchedPatterns);
  DECLARE_BRANCH(double, efficiency);
  DECLARE_BRANCH(double, purity);
  DECLARE_BRANCH(unsigned, nPatHits);
  DECLARE_BRANCH(unsigned, nPatPrecision);
  DECLARE_BRANCH(unsigned, nPatTrigger);
  DECLARE_BRANCH(unsigned, nPatPhi);
  DECLARE_BRANCH(unsigned, nPatPhiOnly);
  DECLARE_BRANCH(unsigned, nPatPhiOnlyMatched);
  DECLARE_BRANCH(int, patSector);
  DECLARE_BRANCH(double, patTheta);
  DECLARE_BRANCH(double, patPhi);
  DECLARE_BRANCH(double, meanNormRes2);

  const GeometryContext gctx{GeometryContext::dangerouslyDefaultConstruct()};
  const TestSeedSelector selector{testCfg.seedCfg};
  const GlobalPatternFinder_t finder{
      makeFinderConfig(selector),
      getDefaultLogger(std::format("GPF_{:}", testCfg.name), gpfLogLvl)};

  const GlobalPatternFinder_t::BeamspotInfo& beamSpot{testCfg.beamSpot};

  RandomEngine engine{testCfg.seed};
  const std::size_t nMuons{testCfg.muonRegions.size()};
  std::size_t nFound{0}, nFakes{0}, nTruePhiOnly{0}, nMatchedPhiOnly{0};
  double totGpfTime{0.};
  const TimePoint_t start = std::chrono::system_clock::now();
  ACTS_INFO("Start test " << testCfg.name << ".");

  for (std::size_t evt = 0; evt < nEvents; ++evt) {
    event = evt;
    std::vector<MuonTruth> muons{};
    Container_t spacePoints{};
    std::vector<TestHit> etaHits{}, phiOnlyHits{};
    std::vector<unsigned> nEtaHitsPerMuon(nMuons, 0),
        nPhiHitsPerMuon(nMuons, 0);
    std::vector<std::array<bool, s_nStations>> stationsPerMuon(nMuons);
    nBkgEtaHits = 0;
    for (const auto [iMuon, region] : enumerate(testCfg.muonRegions)) {
      const MuonTruth& muon{
          muons.emplace_back(generateMuon(engine, testCfg.genCfg, region))};
      const StationCrossings crossings{muonCrossings(muon, testCfg.detCfg)};
      for (const TestHit& hit :
           spawnHits(muon, crossings, static_cast<int>(iMuon), testCfg.detCfg,
                     engine, spacePoints)) {
        if (hit.sp->measuresLoc1()) {
          etaHits.push_back(hit);
          ++nEtaHitsPerMuon[iMuon];
          stationsPerMuon[iMuon][toUnderlying(hit.station)] = true;
        } else {
          phiOnlyHits.push_back(hit);
          ++nPhiHitsPerMuon[iMuon];
        }
      }
      ACTS_DEBUG("Event " << evt << ", muon " << iMuon << ": E = "
                          << muon.energy / 1._GeV << " GeV, eta = " << muon.eta
                          << ", phi = " << muon.phi / 1._degree
                          << ", q = " << muon.charge << " -> "
                          << nEtaHitsPerMuon[iMuon] << " eta hits.");
      if (!testCfg.detCfg.genCfg.addBackground) {
        continue;
      }
      for (const TestHit& hit : spawnBackground(muon, crossings, testCfg.detCfg,
                                                engine, spacePoints)) {
        if (hit.sp->measuresLoc1()) {
          etaHits.push_back(hit);
          ++nBkgEtaHits;
        } else {
          phiOnlyHits.push_back(hit);
        }
      }
    }

    std::vector<GlobalPatternFinder_t::OutputPattern> patterns{};
    if (!etaHits.empty()) {
      const SearchTree_t tree{buildSearchTree(etaHits)};
      const TestPhiProvider phiProvider{&phiOnlyHits};
      const TimePoint_t gpfStart = std::chrono::system_clock::now();
      patterns =
          finder.findPatterns(gctx, tree, selector, phiProvider, beamSpot);
      gpfTime = std::chrono::duration<double, std::micro>(
                    std::chrono::system_clock::now() - gpfStart)
                    .count();
      totGpfTime += gpfTime;
    } else {
      gpfTime = 0.;
    }
    nPatterns = patterns.size();

    /// Count the hits of each muon in each pattern. A pattern is matched to
    /// the muon contributing at least half of its hits.
    std::vector<std::vector<unsigned>> nMatched(
        patterns.size(), std::vector<unsigned>(nMuons, 0));
    std::vector<unsigned> nHits(patterns.size(), 0);
    std::vector<int> patternMuon(patterns.size(), -1);
    for (const auto [iPat, pattern] : enumerate(patterns)) {
      for (const auto& groupHits : pattern.hitsPerGroup) {
        for (const TestHit* hit : groupHits) {
          ++nHits[iPat];
          if (hit->muonIdx >= 0) {
            ++nMatched[iPat][hit->muonIdx];
          }
        }
      }
      for (std::size_t iMuon = 0; iMuon < nMuons; ++iMuon) {
        if (nMatched[iPat][iMuon] >=
            testCfg.patternMatchFraction * nHits[iPat]) {
          patternMuon[iPat] = static_cast<int>(iMuon);
        }
      }
    }
    nFakePatterns = std::ranges::count(patternMuon, -1);
    nFakes += nFakePatterns;

    for (const auto [iMuon, muon] : enumerate(muons)) {
      muonIdx = static_cast<int>(iMuon);
      isBarrelMuon = testCfg.muonRegions[iMuon] == Region::Barrel;
      trueE = muon.energy / 1._GeV;
      trueEta = muon.eta;
      truePhi = muon.phi / 1._degree;
      trueTheta = muon.theta() / 1._degree;
      trueQ = muon.charge;
      trueSector = TestSector{muon.phi}.sector();
      nTrueEtaHits = nEtaHitsPerMuon[iMuon];
      nTruePhiOnlyHits = nPhiHitsPerMuon[iMuon];
      nTrueStations = std::ranges::count(stationsPerMuon[iMuon], true);
      nMatchedPatterns = std::ranges::count(patternMuon, muonIdx);

      /// Select the pattern with the most hits of the muon
      std::optional<std::size_t> bestPattern{};
      for (std::size_t iPat = 0; iPat < patterns.size(); ++iPat) {
        if (nMatched[iPat][iMuon] > 0 &&
            (!bestPattern ||
             nMatched[iPat][iMuon] > nMatched[*bestPattern][iMuon])) {
          bestPattern = iPat;
        }
      }
      const unsigned bestMatched{bestPattern ? nMatched[*bestPattern][iMuon]
                                             : 0u};
      nPatHits = bestPattern ? nHits[*bestPattern] : 0u;
      efficiency = nTrueEtaHits > 0 ? 1. * bestMatched / nTrueEtaHits : 0.;
      purity = nPatHits > 0 ? 1. * bestMatched / nPatHits : 0.;
      found = efficiency >= testCfg.minFoundEfficiency &&
              purity >= testCfg.minFoundPurity;
      nFound += found;
      if (bestPattern) {
        const auto& pattern{patterns[*bestPattern]};
        nPatPrecision = pattern.nPrecisionLayers;
        nPatTrigger = pattern.nTriggerLayers;
        nPatPhi = pattern.nPhiLayers;
        nPatPhiOnly = pattern.phiOnlyHits.size();
        nPatPhiOnlyMatched = std::ranges::count_if(
            pattern.phiOnlyHits,
            [&](const TestHit& hit) { return hit.muonIdx == muonIdx; });
        patSector = pattern.sector.sector();
        patTheta = pattern.patTheta / 1._degree;
        patPhi = pattern.patPhi / 1._degree;
        meanNormRes2 = pattern.meanNormResidual2;
      } else {
        nPatPrecision = nPatTrigger = nPatPhi = nPatPhiOnly = 0;
        nPatPhiOnlyMatched = 0;
        patSector = -1;
        patTheta = patPhi = meanNormRes2 = -999.;
      }
      if (!found) {
        ACTS_DEBUG("Event " << evt << ": muon " << iMuon << " with eta "
                            << trueEta << " and " << nTrueEtaHits
                            << " eta hits not found. Efficiency: " << efficiency
                            << ", purity: " << purity);
      }
      nTruePhiOnly += nTruePhiOnlyHits;
      nMatchedPhiOnly += nPatPhiOnlyMatched;
      outTree->Fill();
    }
  }

  const TimePoint_t end = std::chrono::system_clock::now();
  TestResult result{};
  result.time =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();
  result.foundFraction = static_cast<double>(nFound) / (nEvents * nMuons);
  result.nFakes = nFakes;
  result.gpfTime = totGpfTime / nEvents;
  std::unique_lock guard{writeMutex};
  outFile.WriteObject(outTree.get(), outTree->GetName());
  ACTS_INFO("Test " << testCfg.name << " finished. Found " << nFound << "/"
                    << nEvents * nMuons << " muons, " << nFakes
                    << " patterns not matched to any muon. Attached "
                    << nMatchedPhiOnly << "/" << nTruePhiOnly
                    << " phi-only hits of the muons. It took " << result.time
                    << " ms, pattern finding: " << result.gpfTime
                    << " us/event.");
  return result;
}
#undef DECLARE_BRANCH

/// @brief Configuration of the launched tests and of their checks
struct SuiteConfig {
  /// @brief Number of muons per event. In the Mixed test, barrel and endcap
  ///         muons alternate.
  std::size_t nMuonsPerEvent{2};
  /// @brief Muon energies, one test per energy and event composition
  std::vector<double> energies{5._GeV,  10._GeV,  20._GeV,
                               50._GeV, 100._GeV, 1000._GeV};
  /// @brief Random seeds of the Barrel, Endcap and Mixed tests
  unsigned barrelSeed{1602};
  unsigned endcapSeed{1503};
  unsigned mixedSeed{1701};
  /// @brief Required fraction of found muons for the tests at or above the
  ///         minimum energy
  double minFoundFraction{0.9};
  double minEnergyForCheck{20._GeV};
};

BOOST_AUTO_TEST_SUITE(GlobalPatternFinderTestSuite)

BOOST_AUTO_TEST_CASE(MuonPatternFinding) {
  using namespace std::chrono_literals;
  /// The trees of the tests running in parallel are created concurrently
  ROOT::EnableThreadSafety();

  auto outFile =
      std::make_unique<TFile>("GlobalPatternFinderTest.root", "RECREATE");

  const SuiteConfig suiteCfg{};

  struct LaunchedTest {
    std::string name{};
    double energy{0.};
    std::future<TestResult> result{};
  };
  std::vector<LaunchedTest> tests{};

  auto sendSleep = [&tests]() {
    do {
      std::this_thread::sleep_for(10ms);
    } while (std::ranges::count_if(tests, [](const LaunchedTest& test) {
               return test.result.wait_for(0ms) != std::future_status::ready;
             }) >= nThreads);
  };
  /// @brief Helper lambda to launch the test for a given muon composition
  ///        of the events over the energy points. Each energy point runs in
  ///        its own thread.
  /// @param testName: Name of the test
  /// @param muonRegions: Region of each muon of the event
  /// @param seed: Seed number for the random number generator
  auto launchTest = [&](const std::string& testName,
                        const std::vector<Region>& muonRegions,
                        const unsigned seed) {
    for (const double energy : suiteCfg.energies) {
      TestConfig testCfg{};
      testCfg.name =
          std::format("{:}{:}GeV", testName, static_cast<int>(energy / 1._GeV));
      testCfg.muonRegions = muonRegions;
      testCfg.genCfg.energy = energy;
      testCfg.seed = seed + static_cast<unsigned>(energy / 1._GeV);
      sendSleep();
      tests.emplace_back(testCfg.name, energy,
                         std::async(std::launch::async, [testCfg, &outFile]() {
                           return runPatternTest(testCfg, *outFile);
                         }));
    }
  };
  std::vector<Region> mixedRegions(suiteCfg.nMuonsPerEvent, Region::Barrel);
  for (std::size_t iMuon = 1; iMuon < suiteCfg.nMuonsPerEvent; iMuon += 2) {
    mixedRegions[iMuon] = Region::Endcap;
  }
  launchTest("Barrel",
             std::vector<Region>(suiteCfg.nMuonsPerEvent, Region::Barrel),
             suiteCfg.barrelSeed);
  launchTest("Endcap",
             std::vector<Region>(suiteCfg.nMuonsPerEvent, Region::Endcap),
             suiteCfg.endcapSeed);
  launchTest("Mixed", mixedRegions, suiteCfg.mixedSeed);

  /// Wait until all the tests are completed. The checks are done in the main
  /// thread, since Boost.Test assertions are not thread safe.
  auto timeHisto =
      std::make_unique<TH1D>("TestTimings", "timings; ; time per event [ms]",
                             tests.size(), 0, tests.size());
  auto gpfTimeHisto = std::make_unique<TH1D>(
      "GPFTimings", "pattern finding timings; ; time per event [us]",
      tests.size(), 0, tests.size());
  int bin{1};
  for (LaunchedTest& test : tests) {
    const TestResult result{test.result.get()};
    if (test.energy >= suiteCfg.minEnergyForCheck) {
      BOOST_TEST_CONTEXT("Test " << test.name) {
        BOOST_CHECK_GE(result.foundFraction, suiteCfg.minFoundFraction);
      }
    }
    timeHisto->GetXaxis()->SetBinLabel(bin, test.name.c_str());
    timeHisto->SetBinContent(
        bin, static_cast<double>(result.time) / static_cast<double>(nEvents));
    gpfTimeHisto->GetXaxis()->SetBinLabel(bin, test.name.c_str());
    gpfTimeHisto->SetBinContent(bin, result.gpfTime);
    ++bin;
  }
  outFile->WriteObject(timeHisto.get(), timeHisto->GetName());
  outFile->WriteObject(gpfTimeHisto.get(), gpfTimeHisto->GetName());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

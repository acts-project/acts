// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/ParticleData.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Seeding/GlobalPatternFinder.hpp"
#include "Acts/Utilities/AngleHelpers.hpp"
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
using TimePoint_t = std::chrono::system_clock::time_point;
///////////////////////////////////////////////////////////////////////////////
/// Toy detector description
///////////////////////////////////////////////////////////////////////////////
/// @brief Muon stations. Each station is one group of the pattern topology.
///  There is no transition region: barrel and endcap muons cross only
///  the barrel and only the endcap stations, respectively.
constexpr std::uint8_t s_nStations = 6u;
enum class Station : std::uint8_t { BI = 0, BM, BO, EI, EM, EO };
/// @brief Barrel or endcap region
enum class Region : std::uint8_t { Barrel, Endcap };
/// @brief Station layer
enum class StationLayer : std::uint8_t { Inner, Middle, Outer };
/// @brief Distance of the first tube layer of each station from the beam
///  axis (barrel) or from the transverse plane through the origin
///  (endcap). Each station is about 2.3 m thick.
constexpr std::array<double, s_nStations> s_firstTubeDist{
    4.65_m, 7.0_m, 9.3_m, 7.0_m, 13.0_m, 20.0_m};
/// @brief Station extension along |z| (barrel) or along rho (endcap)
constexpr std::array<double, s_nStations> s_extMin{0._m, 0._m, 0._m,
                                                   1._m, 1._m, 1._m};
constexpr std::array<double, s_nStations> s_extMax{6.5_m, 9.5_m, 12.5_m,
                                                   5._m,  12._m, 12.5_m};
/// @brief Barrel toroid: rhoMin < rho < rhoMax and |z| < endcap zMin
constexpr double s_barrelToroidRhoMin{5.1_m};
constexpr double s_barrelToroidRhoMax{10.5_m};
constexpr double s_barrelToroidField{0.5_T};
/// @brief Endcap toroids: zMin < |z| < zMax and rho < barrel rhoMax
constexpr double s_endcapToroidZMin{9.2_m};
constexpr double s_endcapToroidZMax{12.4_m};
constexpr double s_endcapToroidField{0.5_T};
/// @brief Detector envelope extension
constexpr double s_envelopeMaxRho{14._m};
constexpr double s_envelopeMaxAbsZ{24._m};

/// @brief Whether the station is in the barrel or endcap region
constexpr Region toRegion(const Station st) {
  return (st == Station::BI || st == Station::BM || st == Station::BO)
             ? Region::Barrel
             : Region::Endcap;
}
/// @brief Return the station layer (Inner, Middle, Outer) of a station
constexpr StationLayer toStationLayer(const Station st) {
  return static_cast<StationLayer>(toUnderlying(st) % 3);
}
/// @brief Field strength (along phi) at a given position
double field(const double rho, const double z) {
  const double absZ{std::abs(z)};
  if (rho > s_envelopeMaxRho) {
    return 0.;
  }
  if (absZ >= s_endcapToroidZMin) {
    return absZ < s_endcapToroidZMax ? s_endcapToroidField : 0.;
  }
  return rho > s_barrelToroidRhoMin ? s_barrelToroidField : 0.;
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
/// Mocked GPF ingredients - global pattern finder hit, topology, seed selector
/// and phi-only hit provider
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
  /// @brief Station of the hit
  Station station{Station::BI};
  /// @brief Sector of the hit
  TestSector::Index_t sector{0u};
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
  static constexpr std::uint8_t nGroups{s_nStations};
  /// @brief Compare two hits by measurement layer.
  static std::weak_ordering compareLayers(const TestHit& hit1,
                                          const TestHit& hit2) {
    if (hit1.station != hit2.station) {
      if (toStationLayer(hit1.station) != toStationLayer(hit2.station)) {
        return toStationLayer(hit1.station) <=> toStationLayer(hit2.station);
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
  /// @brief Return the group index of a hit.
  static std::uint8_t groupIndex(const TestHit& hit) {
    return toUnderlying(hit.station);
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
  using Window_t = std::array<double, TestTopology::nGroups>;
  /// @brief Default constructor
  TestSeedSelector() = default;
  /// @brief Whether the hit is good for seeding a pattern.
  bool goodForSeeding(const TestHit& hit) const {
    return toStationLayer(hit.station) != StationLayer::Inner;
  }
  /// @brief Return the theta search window(s) for a given hit.
  const Window_t& thetaSearchWindow(const TestHit& /*hit*/) const {
    return m_window;
  }
 private:
  double m_thetaWindow{3._degree};
  Window_t m_window{filledArray<double, TestTopology::nGroups>(m_thetaWindow)};
};
static_assert(PatternSeedSelector<TestSeedSelector, TestHit, TestTopology>);

/// @brief Provides the phi-only hits from the chambers (station & sector)
///         contributing eta hits to the pattern
struct TestPhiProvider {
  using PatternState_t = GlobalPatternFinder_t::PatternStateAux;
  using PhiHitsPerGroup =
      std::array<std::vector<TestHit>, TestTopology::nGroups>;

  /// @brief Constructor from the phi-only hits of the event
  explicit TestPhiProvider(std::vector<TestHit>&& phiHits)
      : phiOnlyHits{std::move(phiHits)} {}

  /// @brief Return the phi-only hits for each group in the pattern.
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
      for (const TestHit& phiHit : phiOnlyHits) {
        if (TestTopology::groupIndex(phiHit) == group &&
            std::ranges::find(sectors, phiHit.sector) != sectors.end()) {
          result[group].push_back(phiHit);
        }
      }
    }
    return result;
  }
  std::vector<TestHit> phiOnlyHits{};
};
static_assert(
    OnlyPhiHitsProvider<TestPhiProvider, TestHit, TestTopology, TestSector>);

///////////////////////////////////////////////////////////////////////////////
/// Muon generator
///////////////////////////////////////////////////////////////////////////////

/// @brief Configuration of the hit generation in each chamber: the
///  MeasurementGenerator configuration, extended with the options of
///  this test.
using StrawGenCfg_t = ActsTests::MeasurementGenerator::Config;
struct HitGeneratorConfig : StrawGenCfg_t {
  /// @brief MDT + combined trigger strips (RPC/TGC like): 2 gas gaps below
  ///  and above the tubes
  HitGeneratorConfig() {
    createStraws = true;
    twinStraw = false;
    createStrips = true;
    combineSpacePoints = true;
    smearStrips = false;
    discretizeStrips = true;
    stripDirLoc0 = std::vector<Vector3>(2, Vector3::UnitY());
    stripDirLoc1 = std::vector<Vector3>(2, Vector3::UnitX());
  }
  /// @brief Number of phi-only strip layers on each side of the tubes (phi
  ///  strips without an eta partner), outside the combined strip layers
  std::size_t nPhiOnlyLayers{1};
  /// @brief Background: in a few stations crossed by each muon, add an
  ///  independent straight line crossing the muon trajectory.
  std::size_t nBkgStationsPerMuon{2};
  /// @brief Minimum angle between a background line and the
  ///  muon in the bending plane
  double bkgMinAngle{5._degree};
  /// @brief Maximum angle between a background line and the
  ///  muon in the bending plane
  double bkgMaxAngle{40._degree};
  /// @brief Maximum distance of a background line from the muon in the
  ///  crossed station, measured at the station center.
  double bkgMaxOffset{30._cm};
  /// @brief Maximum inclination of a background line and the muon in
  ///  the non-bending plane. Keep it small to stress eta-only separation.
  double bkgMaxTilt{5._degree};
};

/// @brief Truth information of the generated muon
struct MuonTruth {
  double energy{0.};
  double p{0.};
  double eta{0.};
  double phi{0.};
  int charge{1};
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
  const double muonMass{*findMass(PdgParticle::eMuon)};
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
  /// @brief Step size of the trajectory stepping
  static constexpr double s_stepSize{5._mm};
  /// @brief Maximum path length of the trajectory stepping
  static constexpr double s_envelopeMaxPath{40._m};
  explicit ToroidTrajectory(const MuonTruth& muon) {
    /// The angle w.r.t. the transverse plane is pi/2 - theta
    Point point{0., 0., 90._degree - AngleHelpers::thetaFromEta(muon.eta)};
    m_points.push_back(point);
    for (double s = 0.; s < s_envelopeMaxPath; s += s_stepSize) {
      const double kappa{muon.charge * field(point.rho, point.z) / muon.p};
      const double alpha1{point.alpha + kappa * s_stepSize};
      if (std::abs(kappa) < s_epsilon) {
        point.rho += s_stepSize * std::cos(point.alpha);
        point.z += s_stepSize * std::sin(point.alpha);
      } else {
        point.rho += (std::sin(alpha1) - std::sin(point.alpha)) / kappa;
        point.z -= (std::cos(alpha1) - std::cos(point.alpha)) / kappa;
      }
      point.alpha = alpha1;
      /// Stop if the muon leaves the spectrometer or curls back
      if (point.rho < 0. || point.rho > s_envelopeMaxRho ||
          std::abs(point.z) > s_envelopeMaxAbsZ || std::cos(point.alpha) < 0.) {
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

/// @brief Offset between the local z of the MeasurementGenerator and the
///         distance from the beam axis (barrel) or from the transverse plane
///         through the origin (endcap)
double chamberOffset(const Station st) {
  return s_firstTubeDist[toUnderlying(st)] -
         MeasurementGenerator::chamberDistance;
}

/// @brief Local frame of a chamber of the MeasurementGenerator: x along the
///         tubes, y along the bending direction and z orthogonal to the tube
///         layers, pointing outwards. The station plane at the local
///         z = MeasurementGenerator::chamberDistance is placed at the first
///         tube layer distance. Barrel: y = global z, z = sector radial
///         direction. Endcap: y = sector radial direction, z = +- global z.
Transform3 chamberFrame(const Station st, const TestSector& sector,
                        const double zSign) {
  const double offset{chamberOffset(st)};
  Transform3 localToGlobal{Transform3::Identity()};
  if (toRegion(st) == Region::Barrel) {
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

/// @brief Counters of the checks on the MeasurementGenerator output: the test
///        relies on the strip layout, on the phi-only strip convention and on
///        the tube layer positions of the generator. Each counter counts the
///        violations and must stay at zero.
struct GeneratorChecks {
  /// @brief Number of generated chambers
  std::size_t nChambers{0};
  /// @brief Chambers with a number of combined or phi-only strips different
  ///        from the configured strip layers, or with loc1-only strips
  std::size_t nBadStripCounts{0};
  /// @brief Phi-only strips not measuring loc0 or whose half length is not
  ///        stored as the covariance of the bending coordinate
  std::size_t nBadPhiOnlyStrips{0};
  /// @brief Tube hits outside the tube layers of the chamber
  std::size_t nTubesOutside{0};
};

/// @brief Create the space points of a straight line crossing a station. The
///         MeasurementGenerator is run on the line expressed in the local
///         chamber frame.
///  @param globPos: Point on the line, in global coordinates
///  @param globDir: Direction of the line, in global coordinates
///  @param muonIdx: Index of the generating muon (-1 for background)
///  @param cfg: Configuration of the hit generation
///  @param outSpacePoints: Container owning the generated space points
///  @param outHits: Hits wrapping the generated space points
///  @param checks: Checks on the generator output
void spawnInStation(const Vector3& globPos, const Vector3& globDir,
                    const Station station, const TestSector& sector,
                    const double zSign, const int muonIdx,
                    const HitGeneratorConfig& cfg, RandomEngine& engine,
                    Container_t& outSpacePoints, std::vector<TestHit>& outHits,
                    GeneratorChecks& checks) {
  const auto st{toUnderlying(station)};
  const Transform3 localToGlobal{chamberFrame(station, sector, zSign)};
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
  StrawGenCfg_t genCfg{cfg};
  genCfg.stripDirLoc0.resize(genCfg.stripDirLoc1.size() + cfg.nPhiOnlyLayers,
                             Vector3::UnitY());
  genCfg.localToGlobal = localToGlobal;
  Container_t spacePoints{
      MeasurementGenerator::spawn(line, 0., engine, genCfg, logger())};

  /// Check the generator output, before the station extension cut: one
  /// combined strip per loc1 strip layer and one phi-only strip per extra loc0
  /// strip layer, below and above the tubes, and the tube hits within the
  /// tube layers
  using CovIdx =
      Acts::Experimental::detail::CompSpacePointAuxiliaries::ResidualIdx;
  ++checks.nChambers;
  std::size_t nCombined{0}, nPhiOnly{0}, nLoc1Only{0};
  /// Tolerance on the tube layer positions, which the generator computes with
  /// std::sin instead of the constexpr sqrt(3)
  constexpr double tubeZTolerance{1._um};
  const double lastTubeZ{2. * MeasurementGenerator::tubeLayersCentreZ() -
                         MeasurementGenerator::chamberDistance};
  for (const auto& sp : spacePoints) {
    if (sp->isStraw()) {
      const double z{sp->localPosition().z()};
      if (z < MeasurementGenerator::chamberDistance - tubeZTolerance ||
          z > lastTubeZ + tubeZTolerance) {
        ++checks.nTubesOutside;
      }
    } else if (sp->measuresLoc0() && sp->measuresLoc1()) {
      ++nCombined;
    } else if (sp->measuresLoc0()) {
      ++nPhiOnly;
    } else {
      ++nLoc1Only;
    }
  }
  if (nCombined != 2 * genCfg.stripDirLoc1.size() ||
      nPhiOnly != 2 * cfg.nPhiOnlyLayers || nLoc1Only != 0) {
    ++checks.nBadStripCounts;
  }

  for (auto& sp : spacePoints) {
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
      const double ySign{
          toRegion(station) == Region::Barrel && locTrackPos.y() < 0. ? -1.
                                                                      : 1.};
      sp->updatePosition(Vector3{sp->localPosition().x(),
                                 ySign * 0.5 * (s_extMin[st] + s_extMax[st]),
                                 sp->localPosition().z()});
      const double halfLength{0.5 * (s_extMax[st] - s_extMin[st])};
      sp->setStripHalfLength(halfLength);
      if (!sp->measuresLoc0() ||
          sp->covariance()[toUnderlying(CovIdx::bending)] !=
              square(halfLength)) {
        ++checks.nBadPhiOnlyStrips;
      }
    }
    hit.position = localToGlobal * sp->localPosition();
    /// The generated chambers are unbounded: apply the station extension. For
    /// phi-only strips, the track position on the strip plane is checked.
    const Vector3 extPos{isPhiOnly ? Vector3{localToGlobal * locTrackPos}
                                   : hit.position};
    const double ext{toRegion(station) == Region::Barrel
                         ? std::abs(extPos.z())
                         : VectorHelpers::perp(extPos)};
    if (ext < s_extMin[st] || ext > s_extMax[st]) {
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
    hit.layerDistance = chamberOffset(station) + sp->localPosition().z();
    hit.muonIdx = muonIdx;
    outHits.push_back(hit);
    outSpacePoints.push_back(std::move(sp));
  }
}

/// @brief Distance of the centre plane of a station from the beam axis
///         (barrel) or from the transverse plane through the origin (endcap)
double stationCentreDist(const Station station) {
  return chamberOffset(station) + MeasurementGenerator::tubeLayersCentreZ();
}

/// @brief Point and tangent of the muon trajectory at the centre plane of each
///         station of its sector (nullopt if the station is not crossed)
using StationCrossings =
    std::array<std::optional<std::pair<Vector3, Vector3>>, s_nStations>;
StationCrossings muonCrossings(const MuonTruth& muon) {
  const ToroidTrajectory trajectory{muon};
  const Vector3 muonRadialDir{std::cos(muon.phi), std::sin(muon.phi), 0.};
  const TestSector sector{muon.phi};
  const double cosDPhi{std::cos(sector.deltaPhi(muon.phi))};
  StationCrossings crossings{};
  for (std::size_t st = 0; st < s_nStations; ++st) {
    const auto station{static_cast<Station>(st)};
    /// The muon moves in the half plane at its own phi, which crosses the
    /// barrel planes of the sector at rho * cos(dPhi) = centreDist.
    const double centreDist{stationCentreDist(station)};
    const auto point{trajectory.crossing([&](const ToroidTrajectory::Point& p) {
      return toRegion(station) == Region::Barrel ? p.rho * cosDPhi - centreDist
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
                               const int muonIdx, const HitGeneratorConfig& cfg,
                               RandomEngine& engine,
                               Container_t& outSpacePoints,
                               GeneratorChecks& checks) {
  const TestSector sector{muon.phi};
  const double zSign{zSide(muon)};
  std::vector<TestHit> hits{};
  for (const auto [st, crossing] : enumerate(crossings)) {
    if (crossing) {
      spawnInStation(crossing->first, crossing->second,
                     static_cast<Station>(st), sector, zSign, muonIdx, cfg,
                     engine, outSpacePoints, hits, checks);
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
                                     const HitGeneratorConfig& cfg,
                                     RandomEngine& engine,
                                     Container_t& outSpacePoints,
                                     GeneratorChecks& checks) {
  const HitGeneratorConfig& genCfg{cfg};
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
                   -1, cfg, engine, outSpacePoints, hits, checks);
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
  /// @brief Configuration of the hit generation in each chamber
  HitGeneratorConfig hitCfg{};
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
  /// @brief Fraction of the phi-only hits of the muons attached to their
  ///        best pattern
  double phiOnlyFraction{0.};
  /// @brief Whether the background was added
  bool withBackground{false};
  /// @brief Checks on the generator output
  GeneratorChecks genChecks{};
};

GlobalPatternFinder_t::Config makeFinderConfig() {
  /// Default values of the MuonGlobalPatternFindingAlg
  GlobalPatternFinder_t::Config cfg{};
  cfg.nResidualSigma = 3.;
  cfg.lowConfidenceResSigma = 50.;
  cfg.maxThetaSizeOverlap = 6._degree;
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
  const TestSeedSelector selector{};
  const GlobalPatternFinder_t finder{
      makeFinderConfig(),
      getDefaultLogger(std::format("GPF_{:}", testCfg.name), gpfLogLvl)};

  const GlobalPatternFinder_t::BeamspotInfo& beamSpot{testCfg.beamSpot};

  RandomEngine engine{testCfg.seed};
  const std::size_t nMuons{testCfg.muonRegions.size()};
  std::size_t nFound{0}, nFakes{0}, nTruePhiOnly{0}, nMatchedPhiOnly{0};
  GeneratorChecks genChecks{};
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
      const StationCrossings crossings{muonCrossings(muon)};
      for (const TestHit& hit :
           spawnHits(muon, crossings, static_cast<int>(iMuon), testCfg.hitCfg,
                     engine, spacePoints, genChecks)) {
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
      if (testCfg.hitCfg.nBkgStationsPerMuon == 0) {
        continue;
      }
      for (const TestHit& hit :
           spawnBackground(muon, crossings, testCfg.hitCfg, engine, spacePoints,
                           genChecks)) {
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
      const TestPhiProvider phiProvider{std::move(phiOnlyHits)};
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
      trueTheta = AngleHelpers::thetaFromEta(muon.eta) / 1._degree;
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
  result.phiOnlyFraction =
      nTruePhiOnly > 0 ? static_cast<double>(nMatchedPhiOnly) / nTruePhiOnly
                       : 0.;
  result.withBackground = testCfg.hitCfg.nBkgStationsPerMuon > 0;
  result.genChecks = genChecks;
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
  /// @brief Required fraction of the muon phi-only hits attached to their
  ///        pattern, for the tests at or above the minimum energy without
  ///        background (with background, the phi-only strips of the background
  ///        lines crossing the muon compete with the muon ones)
  double minPhiOnlyFraction{0.95};
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
    BOOST_TEST_CONTEXT("Test " << test.name) {
      /// The generator output is what the test expects
      BOOST_CHECK_GT(result.genChecks.nChambers, 0u);
      BOOST_CHECK_EQUAL(result.genChecks.nBadStripCounts, 0u);
      BOOST_CHECK_EQUAL(result.genChecks.nBadPhiOnlyStrips, 0u);
      BOOST_CHECK_EQUAL(result.genChecks.nTubesOutside, 0u);
      if (test.energy >= suiteCfg.minEnergyForCheck) {
        BOOST_CHECK_GE(result.foundFraction, suiteCfg.minFoundFraction);
        if (!result.withBackground) {
          BOOST_CHECK_GE(result.phiOnlyFraction, suiteCfg.minPhiOnlyFraction);
        }
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

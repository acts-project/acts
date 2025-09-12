// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/CompositeSpacePointLineFitter.hpp"
#include "Acts/Seeding/CompositeSpacePointLineSeeder.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <format>
#include <random>
#include <ranges>
#include <span>

#include <TFile.h>
#include <TTree.h>

using namespace Acts;
using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;
using namespace Acts::UnitLiterals;
using namespace Acts::detail::LineHelper;
using namespace Acts::PlanarHelper;
using namespace Acts::VectorHelpers;

using RandomEngine = std::mt19937;
using uniform = std::uniform_real_distribution<double>;
using normal_t = std::normal_distribution<double>;
using Line_t = CompSpacePointAuxiliaries::Line_t;
using ResidualIdx = CompSpacePointAuxiliaries::ResidualIdx;
using FitParIndex = CompSpacePointAuxiliaries::FitParIndex;
using ParamVec_t = CompositeSpacePointLineFitter::ParamVec_t;

constexpr auto logLvl = Acts::Logging::Level::INFO;
constexpr std::size_t nEvents = 10000;

ACTS_LOCAL_LOGGER(getDefaultLogger("StrawLineFitterTest", logLvl));

namespace Acts::Test {

class FitTestSpacePoint {
 public:
  /// @brief Constructor for straw wires
  /// @param pos: Position of the wire
  /// @param driftR: Straw drift radius
  /// @param driftRUncert: Uncertainty on the drift radius uncertainty
  /// @param mLoc0: Flag toggling whether the position along the wire is also provided
  FitTestSpacePoint(const Vector3& pos, const double driftR,
                    const double driftRUncert, const bool mLoc0 = false)
      : m_position{pos}, m_driftR{driftR}, m_measLoc0{mLoc0} {
    m_covariance[toUnderlying(ResidualIdx::bending)] =
        Acts::square(driftRUncert);
  }

  /// @brief Constructor for strip measurements
  FitTestSpacePoint(const Vector3& stripPos, const Vector3& stripDir,
                    const Vector3& toNext, const double uncertLoc0,
                    const double uncertLoc1)
      : m_position{stripPos},
        m_sensorDir{stripDir},
        m_toNextSen{toNext},
        m_measLoc0{uncertLoc0 > 0.},
        m_measLoc1{uncertLoc1 > 0.} {
    using enum ResidualIdx;
    m_covariance[toUnderlying(nonBending)] = Acts::square(uncertLoc0);
    m_covariance[toUnderlying(bending)] = Acts::square(uncertLoc1);
  }
  /// @brief Position of the space point
  const Vector3& localPosition() const { return m_position; }
  /// @brief Wire direction
  const Vector3& sensorDirection() const { return m_sensorDir; }
  /// @brief To next sensor in the plane
  const Vector3& toNextSensor() const { return m_toNextSen; }
  /// @brief To next straw layer
  const Vector3& planeNormal() const { return m_planeNorm; }
  /// @brief Measurement's radius
  double driftRadius() const { return m_driftR.value_or(0.); }
  /// @brief Measurement's covariance
  const std::array<double, 3>& covariance() const { return m_covariance; }
  /// @brief Time of record
  double time() const { return m_time.value_or(0.); }
  /// @brief All measurements are straws
  bool isStraw() const { return m_driftR.has_value(); }
  /// @brief Check whether the space point has a time value
  bool hasTime() const { return m_time.has_value(); }
  /// @brief Check whether the space point measures the non-bending direction
  bool measuresLoc0() const { return m_measLoc0; }
  /// @brief Check whether the space point measures the bending direction
  bool measuresLoc1() const { return m_measLoc1 || isStraw(); }
  /// @brief Sets the straw tube's drift radius
  void updateDriftR(const double updatedR) { m_driftR = updatedR; }
  /// @brief Updates the position of the space point
  void updatePosition(const Vector3& newPos) { m_position = newPos; }

 private:
  Vector3 m_position{Vector3::Zero()};
  Vector3 m_sensorDir{Vector3::UnitX()};
  Vector3 m_toNextSen{Vector3::UnitY()};
  Vector3 m_planeNorm{m_sensorDir.cross(m_toNextSen).normalized()};
  std::optional<double> m_driftR{std::nullopt};
  std::optional<double> m_time{std::nullopt};
  std::array<double, 3> m_covariance{filledArray<double, 3>(0)};
  bool m_measLoc0{false};
  bool m_measLoc1{false};
};

static_assert(CompositeSpacePoint<FitTestSpacePoint>);

using Container_t = std::vector<std::shared_ptr<FitTestSpacePoint>>;

static_assert(CompositeSpacePointContainer<Container_t>);
class SpCalibrator {
 public:
  static double driftUncert(const double r) {
    return 0.2_mm / (1._mm + Acts::square(r)) + 0.1_mm;
  }
  /// @brief Calibrate a set of straw measurements using the best known estimate on a straight line track
  /// @param ctx: Calibration context (Needed by conept interface)
  /// @param trackPos: Position of the track at z=0.
  /// @param trackDir: Direction of the track in the local frame
  /// @param timeOffSet: Offset in the time of arrival (To be implemented)
  /// @param uncalibCont: Uncalibrated composite space point container
  Container_t calibrate(const Acts::CalibrationContext& /*ctx*/,
                        const Vector3& trackPos, const Vector3& trackDir,
                        const double /*timeOffSet*/,
                        const Container_t& uncalibCont) const {
    Container_t calibMeas{};
    for (const auto& sp : uncalibCont) {
      if (!sp->measuresLoc0() || !sp->measuresLoc1()) {
        /// Estimate the best position along the sensor
        auto bestPos = lineIntersect(trackPos, trackDir, sp->localPosition(),
                                     sp->sensorDirection());
        sp->updatePosition(bestPos.position());
      }
      if (sp->isStraw()) {
        sp->updateDriftR(Acts::abs(sp->driftRadius()));
      }
    }
    return uncalibCont;
  }
  /// @brief Updates the sign of the Straw's drift radii indicating that they are on the left (-1)
  ///        or right side (+1) of the track line
  void updateSigns(const Vector3& trackPos, const Vector3& trackDir,
                   Container_t& measurements) const {
    auto signs =
        CompSpacePointAuxiliaries::strawSigns(trackPos, trackDir, measurements);
    /// The signs have the same size as the measurement container
    for (std::size_t s = 0; s < signs.size(); ++s) {
      /// Take care to not turn strips into straws by accident
      if (measurements[s]->isStraw()) {
        measurements[s]->updateDriftR(
            Acts::abs(measurements[s]->driftRadius()) * signs[s]);
      }
    }
  }
};
/// Ensure that the Test space point calibrator satisfies the calibrator concept
static_assert(
    CompositeSpacePointCalibrator<SpCalibrator, Container_t, Container_t>);

/// @brief Generates a random straight line
/// @param engine Random number sequence to draw the parameters from
/// @return A Line object instantiated with the generated parameters
Line_t generateLine(RandomEngine& engine) {
  using ParIndex = Line_t::ParIndex;
  Line_t::ParamVector linePars{};
  linePars[toUnderlying(ParIndex::phi)] = 90._degree;
  // linePars[toUnderlying(ParIndex::phi)] =
  //     std::uniform_real_distribution{-120_degree, 120_degree}(engine);
  linePars[toUnderlying(ParIndex::x0)] =
      std::uniform_real_distribution{-5000., 5000.}(engine);
  linePars[toUnderlying(ParIndex::y0)] =
      std::uniform_real_distribution{-5000., 5000.}(engine);
  linePars[toUnderlying(ParIndex::theta)] =
      std::uniform_real_distribution{0.1_degree, 179.9_degree}(engine);
  if (Acts::abs(linePars[toUnderlying(ParIndex::theta)] - 90._degree) <
      0.2_degree) {
    return generateLine(engine);
  }
  Line_t line{linePars};
  if (line.direction().z() < 0.) {
    linePars[toUnderlying(ParIndex::theta)] = theta(-line.direction());
    linePars[toUnderlying(ParIndex::phi)] = phi(-line.direction());
    line.updateParameters(linePars);
  }

  ACTS_DEBUG(
      "\n\n\n"
      << __func__ << "() " << __LINE__ << " - Generated parameters "
      << std::format("theta: {:.2f}, phi: {:.2f}, y0: {:.1f}, x0: {:.1f}",
                     linePars[toUnderlying(ParIndex::theta)] / 1._degree,
                     linePars[toUnderlying(ParIndex::phi)] / 1._degree,
                     linePars[toUnderlying(ParIndex::y0)],
                     linePars[toUnderlying(ParIndex::x0)])
      << " --> " << toString(line.position()) << " + "
      << toString(line.direction()));

  return line;
}

class MeasurementGenerator {
 public:
  struct Config {
    /// @brief Create straw measurements
    bool createStraws{true};
    /// @brief Smear the straw radius
    bool smearRadius{true};
    /// @brief Create strip measurements
    bool createStrips{true};
    /// @brief Discretize strips
    bool smearStrips{true};
    /// @brief Combine the two strip measurements to a single space point
    bool combineSpacePoints{false};
    /// @brief Create strip measurements constraining Loc0
    bool createStripsLoc0{true};
    /// @brief Create strip measurements constraining Loc1
    bool createStripsLoc1{true};
    /// @brief Pitch between two loc0 strips
    double stripPitchLoc0{4._cm};
    /// @brief Pitch between two loc1 strips
    double stripPitchLoc1{3._cm};
    /// @brief Vector
    Vector3 stripDirLoc1{Vector3::UnitX()};
    /// @brief
    Vector3 stripDirLoc0{Vector3::UnitY()};
  };
  /// @brief Extrapolate the straight line track through the straw layers to
  ///        evaluate which tubes were crossed by the track. The straw layers
  ///        are staggered in the z-direction and each layer expands in y. To
  ///        estimate, which tubes were crossed, the track is extrapolated to
  ///        the z-plane below & above the straw wires. Optionally, the true
  ///        drift-radius can be smeared assuming a Gaussian with a drift-radius
  ///        dependent uncertainty.
  /// @param line: The track to extrapolate
  /// @param engine: Random number generator to smear the drift radius
  /// @param smearRadius: If true, the drift radius is smeared with a Gaussian
  /// @param createStrips: If true the strip measurements are created
  static Container_t spawn(const Line_t& line, const double /*t0*/,
                           RandomEngine& engine, const Config& genCfg) {
    /// Direction vector to go a positive step in the tube honeycomb grid
    const Vector3 posStaggering{0., std::cos(60._degree), std::sin(60._degree)};
    /// Direction vector to go a negative step in the tube honeycomb grid
    const Vector3 negStaggering{0., -std::cos(60._degree),
                                std::sin(60._degree)};
    /// Number of tube layers per multilayer
    constexpr std::size_t nLayersPerMl = 8;
    /// Number of overall tubelayers
    constexpr std::size_t nTubeLayers = nLayersPerMl * 2;
    /// Position in z of the first tube layer
    constexpr double chamberDistance = 3._m;
    /// Radius of each straw
    constexpr double tubeRadius = 15._mm;
    /// Distance between the first <nLayersPerMl> layers and the second pack
    constexpr double tubeLayerDist = 1.2_m;
    /// Distance between the tube multilayers and to the first strip layer
    constexpr double tubeStripDist = 30._cm;
    /// Distance between two strip layers
    constexpr double stripLayDist = 0.5_cm;
    /// Number of strip layers on each side
    constexpr std::size_t nStripLay = 4;

    std::array<Vector3, nTubeLayers> tubePositions{
        filledArray<Vector3, nTubeLayers>(chamberDistance * Vector3::UnitZ())};
    /// Fill the positions of the reference tubes 1
    for (std::size_t l = 1; l < nTubeLayers; ++l) {
      const Vector3& layStag{l % 2 == 1 ? posStaggering : negStaggering};
      tubePositions[l] = tubePositions[l - 1] + 2. * tubeRadius * layStag;

      if (l == nLayersPerMl) {
        tubePositions[l] += tubeLayerDist * Vector3::UnitZ();
      }
    }
    /// Print the staggering
    ACTS_DEBUG("##############################################");

    for (std::size_t l = 0; l < nTubeLayers; ++l) {
      ACTS_DEBUG("  *** " << (l + 1) << " - " << toString(tubePositions[l]));
    }
    ACTS_DEBUG("##############################################");

    Container_t measurements{};
    /// Extrapolate the track to the z-planes of the tubes and determine which
    /// tubes were actually hit
    if (genCfg.createStraws) {
      for (const auto& stag : tubePositions) {
        auto planeExtpLow =
            intersectPlane(line.position(), line.direction(), Vector3::UnitZ(),
                           stag.z() - tubeRadius);
        auto planeExtpHigh =
            intersectPlane(line.position(), line.direction(), Vector3::UnitZ(),
                           stag.z() + tubeRadius);
        ACTS_DEBUG("Extrapolated to plane "
                   << toString(planeExtpLow.position()) << " "
                   << toString(planeExtpHigh.position()));

        const auto dToFirstLow = static_cast<int>(std::ceil(
            (planeExtpLow.position().y() - stag.y()) / (2. * tubeRadius)));
        const auto dToFirstHigh = static_cast<int>(std::ceil(
            (planeExtpHigh.position().y() - stag.y()) / (2. * tubeRadius)));
        /// Does the track go from left to right or right to left?
        const int dT = dToFirstHigh > dToFirstLow ? 1 : -1;
        /// Loop over the candidate tubes and check each one whether the track
        /// actually crossed them. Then generate the circle and optionally smear
        /// the radius
        for (int tN = dToFirstLow - dT; tN != dToFirstHigh + 2 * dT; tN += dT) {
          const Vector3 tube = stag + 2. * tN * tubeRadius * Vector3::UnitY();
          const double rad = Acts::abs(Acts::detail::LineHelper::signedDistance(
              tube, Vector3::UnitX(), line.position(), line.direction()));
          if (rad > tubeRadius) {
            continue;
          }
          ACTS_DEBUG("Tube position: " << toString(tube)
                                       << ", radius: " << rad);

          const double smearedR =
              genCfg.smearRadius
                  ? Acts::abs(
                        normal_t{rad, SpCalibrator::driftUncert(rad)}(engine))
                  : rad;
          if (smearedR > tubeRadius) {
            continue;
          }
          measurements.emplace_back(std::make_unique<FitTestSpacePoint>(
              tube, smearedR, SpCalibrator::driftUncert(smearedR)));
        }
      }
    }
    if (genCfg.createStrips) {
      ///@brief Helper function to discretize the extrapolated position on the strip plane to actual strips.
      ///       The function returns the local coordinate in the picked plane
      ///       projection
      /// @param extPos: Extrapolated position of the straight line in the plane
      /// @param loc1: Boolean to fetch either the loc1 projection or the loc0 projection
      auto discretize = [&genCfg](const Vector3& extPos, const bool loc1) {
        const double pitch =
            loc1 ? genCfg.stripPitchLoc1 : genCfg.stripPitchLoc0;
        assert(pitch > 0.);

        const Vector3& stripDir =
            loc1 ? genCfg.stripDirLoc1 : genCfg.stripDirLoc0;
        const Vector3 stripNormal = stripDir.cross(Vector3::UnitZ());
        const double dist = stripNormal.dot(extPos);
        return genCfg.smearStrips
                   ? pitch * static_cast<int>(std::ceil(dist / pitch))
                   : dist;
      };
      /// Calculate the strip measurement's covariance
      const double stripCovLoc0 =
          Acts::square(genCfg.stripPitchLoc0) / std::sqrt(12.);
      const double stripCovLoc1 =
          Acts::square(genCfg.stripPitchLoc1) / std::sqrt(12.);

      for (std::size_t sL = 0; sL < nStripLay; ++sL) {
        const double distInZ = tubeStripDist + sL * stripLayDist;
        const double planeLow = tubePositions.front().z() - distInZ;
        const double planeHigh = tubePositions.back().z() + distInZ;

        for (const double plane : {planeLow, planeHigh}) {
          const auto extp = intersectPlane(line.position(), line.direction(),
                                           Vector3::UnitZ(), plane);

          if (genCfg.combineSpacePoints) {
            const Vector3 extpPos{discretize(extp.position(), false),
                                  discretize(extp.position(), true), plane};
            measurements.emplace_back(std::make_unique<FitTestSpacePoint>(
                extpPos, genCfg.stripDirLoc0, genCfg.stripDirLoc1, stripCovLoc0,
                stripCovLoc1));
          } else {
            if (genCfg.createStripsLoc0) {
              const Vector3 extpPos{discretize(extp.position(), false), 0.,
                                    plane};
              measurements.emplace_back(std::make_unique<FitTestSpacePoint>(
                  extpPos, genCfg.stripDirLoc0,
                  genCfg.stripDirLoc0.cross(Vector3::UnitZ()), stripCovLoc0,
                  0.));
            }
            if (genCfg.createStripsLoc1) {
              const Vector3 extpPos{0., discretize(extp.position(), true),
                                    plane};
              measurements.emplace_back(std::make_unique<FitTestSpacePoint>(
                  extpPos, genCfg.stripDirLoc1.cross(Vector3::UnitZ()),
                  genCfg.stripDirLoc1, 0., stripCovLoc1));
            }
          }
        }
      }
    }
    ACTS_DEBUG("Track hit in total " << measurements.size() << " sensors.");
    std::ranges::sort(measurements, [&line](const auto& spA, const auto& spB) {
      return line.direction().dot(spA->localPosition() - line.position()) <
             line.direction().dot(spB->localPosition() - line.position());
    });
    return measurements;
  }
};
/// @brief Construct the start parameters from the hit container && the true trajectory
/// @param line: True trajectory to pick-up the correct left/right ambiguity for the straw seed hits
/// @param hits: List of measurements to be used for fitting
ParamVec_t startParameters(const Line_t& line, const Container_t& hits) {
  ParamVec_t pars{};

  double tanPhi{0.};
  double tanTheta{0.};
  /// Setup the seed parameters in x0 && phi
  auto firstPhi = std::ranges::find_if(hits, [](const auto& sp) {
    return !sp->isStraw() && sp->measuresLoc0();
  });
  auto lastPhi = std::ranges::find_if(
      std::ranges::reverse_view(hits),
      [](const auto& sp) { return !sp->isStraw() && sp->measuresLoc0(); });

  if (firstPhi != hits.end() && lastPhi != hits.rend()) {
    const Vector3 firstToLastPhi =
        (**lastPhi).localPosition() - (**firstPhi).localPosition();
    tanPhi = firstToLastPhi.x() / firstToLastPhi.z();
    /// -> x = tanPhi * z + x_{0} ->
    pars[toUnderlying(FitParIndex::x0)] =
        (**lastPhi).localPosition().x() -
        (**lastPhi).localPosition().z() * tanPhi;
  }
  /// Setup the seed parameters in y0 && theta
  auto firstTube =
      std::ranges::find_if(hits, [](const auto& sp) { return sp->isStraw(); });
  auto lastTube =
      std::ranges::find_if(std::ranges::reverse_view(hits),
                           [](const auto& sp) { return sp->isStraw(); });

  if (firstTube != hits.end() && lastTube != hits.rend()) {
    const int signFirst =
        CompSpacePointAuxiliaries::strawSign(line, **firstTube);
    const int signLast = CompSpacePointAuxiliaries::strawSign(line, **lastTube);

    auto seedPars = CompositeSpacePointLineSeeder::constructTangentLine(
        **lastTube, **firstTube,
        CompositeSpacePointLineSeeder::encodeAmbiguity(signLast, signFirst));
    tanTheta = std::tan(seedPars.theta);
    pars[toUnderlying(FitParIndex::y0)] = seedPars.y0;
  }

  const Vector3 seedDir = makeDirectionFromAxisTangents(tanPhi, tanTheta);
  pars[toUnderlying(FitParIndex::theta)] = theta(seedDir);
  pars[toUnderlying(FitParIndex::phi)] = phi(seedDir);
  return pars;
}

BOOST_AUTO_TEST_SUITE(FastStrawLineFitTests)

BOOST_AUTO_TEST_CASE(SeedTangents) {
  RandomEngine engine{1602};
  constexpr double tolerance = 1.e-3;
  return;

  using Seeder = CompositeSpacePointLineSeeder;
  using SeedAux = CompSpacePointAuxiliaries;
  using enum Seeder::TangentAmbi;

  for (std::size_t evt = 0; evt < nEvents; ++evt) {
    const auto line = generateLine(engine);
    auto testTubes = generateMeasurements(line, 0._ns, engine, false, false);
    const double lineTanTheta = line.direction().y() / line.direction().z();
    const double lineY0 = line.position().y();
    for (std::size_t m1 = testTubes.size() - 1; m1 > 0; --m1) {
      for (std::size_t m2 = 0; m2 < m1; ++m2) {
        const auto& bottomTube = *testTubes[m2];
        const auto& topTube = *testTubes[m1];
        const int signTop = SeedAux::strawSign(line, topTube);
        const int signBot = SeedAux::strawSign(line, bottomTube);
        const auto trueAmbi = Seeder::encodeAmbiguity(signTop, signBot);

        ACTS_DEBUG(__func__ << "() " << __LINE__ << " - "
                            << std::format("bottom tube @ {:}, r: {:.3f}({:})",
                                           toString(bottomTube.localPosition()),
                                           (signBot * bottomTube.driftRadius()),
                                           signBot > 0 ? "R" : "L")
                            << ", "
                            << std::format("top tube @ {:}, r: {:.3f} ({:})",
                                           toString(topTube.localPosition()),
                                           (signTop * topTube.driftRadius()),
                                           signTop > 0 ? "R" : "L"));

        bool seenTruePars{false};
        for (const auto ambi : {LL, RL, LR, RR}) {
          const auto pars =
              Seeder::constructTangentLine(topTube, bottomTube, ambi);

          const bool isTruePars =
              Acts::abs(std::tan(pars.theta) - lineTanTheta) < tolerance &&
              Acts::abs(lineY0 - pars.y0) < tolerance;
          seenTruePars |= isTruePars;
          ACTS_VERBOSE(__func__
                       << "() " << __LINE__ << " - Test ambiguity "
                       << CompositeSpacePointLineSeeder::toString(ambi)
                       << " -> "
                       << std::format("theta: {:.3f}, ", pars.theta / 1._degree)
                       << std::format("tanTheta: {:.3f}, ",
                                      std::tan(pars.theta))
                       << std::format("y0: {:.3f}", pars.y0)
                       << (!isTruePars ? (ambi == trueAmbi ? " xxxxxx" : "")
                                       : " <-------"));
          const Vector3 seedPos = pars.y0 * Vector3::UnitY();
          const Vector3 seedDir =
              makeDirectionFromAxisTangents(0., std::tan(pars.theta));

          const double chi2Top = SeedAux::chi2Term(seedPos, seedDir, topTube);
          const double chi2Bot =
              SeedAux::chi2Term(seedPos, seedDir, bottomTube);
          ACTS_VERBOSE(__func__ << "() " << __LINE__
                                << " - Resulting chi2 top: " << chi2Top
                                << ", bottom: " << chi2Bot);
          BOOST_CHECK_LE(chi2Top, 1.e-17);
          BOOST_CHECK_LE(chi2Bot, 1.e-17);
        }
        BOOST_CHECK_EQUAL(seenTruePars, true);

        const auto seedPars =
            Seeder::constructTangentLine(topTube, bottomTube, trueAmbi);
        /// Construct line parameters
        ACTS_DEBUG(__func__
                   << "() " << __LINE__ << " - Line tan theta: " << lineTanTheta
                   << ", reconstructed theta: " << std::tan(seedPars.theta)
                   << ", line y0: " << lineY0
                   << ", reconstructed y0: " << seedPars.y0);
        BOOST_CHECK_CLOSE(std::tan(seedPars.theta), lineTanTheta, tolerance);
        BOOST_CHECK_CLOSE(seedPars.y0, lineY0, tolerance);
      }
    }
  }
}

#define DECLARE_BRANCH(dType, bName) \
  dType bName{};                     \
  outTree->Branch(#bName, &bName);

BOOST_AUTO_TEST_CASE(SimpleLineFit) {
  RandomEngine engine{1503};

  using Fitter_t = CompositeSpacePointLineFitter;
  using FitOpts_t = Fitter_t::FitOptions<Container_t, SpCalibrator>;

  Fitter_t::Config cfg{};
  cfg.useHessian = false;

  Fitter_t fitter{cfg, getDefaultLogger("LineFitter", logLvl)};

  auto outFile =
      std::make_unique<TFile>("StrawLineFitterTest.root", "RECREATE");
  auto outTree = std::make_unique<TTree>("SimpleLineFit", "MonitorTree");

  DECLARE_BRANCH(double, trueY0);
  DECLARE_BRANCH(double, trueX0);
  DECLARE_BRANCH(double, trueTheta);
  DECLARE_BRANCH(double, truePhi);

  DECLARE_BRANCH(double, recoY0);
  DECLARE_BRANCH(double, recoX0);
  DECLARE_BRANCH(double, recoTheta);
  DECLARE_BRANCH(double, recoPhi);

  DECLARE_BRANCH(double, covY0);
  DECLARE_BRANCH(double, covX0);
  DECLARE_BRANCH(double, covTheta);
  DECLARE_BRANCH(double, covPhi);

  DECLARE_BRANCH(double, chi2);
  DECLARE_BRANCH(unsigned, nIter);
  DECLARE_BRANCH(unsigned, nDoF);

  auto fillPars = [](const auto pars, double& y0, double& x0, double& theta,
                     double& phi) {
    y0 = pars[toUnderlying(FitParIndex::y0)];
    x0 = pars[toUnderlying(FitParIndex::x0)];
    theta = pars[toUnderlying(FitParIndex::theta)];
    phi = pars[toUnderlying(FitParIndex::phi)];
  };

  auto calibrator = std::make_unique<SpCalibrator>();
  for (std::size_t evt = 0; evt < nEvents; ++evt) {
    const auto line = generateLine(engine);
    fillPars(line.parameters(), trueY0, trueX0, trueTheta, truePhi);
    const double t0 = uniform{-50._ns, 50._ns}(engine);

    FitOpts_t fitOpts{};
    fitOpts.calibrator = calibrator.get();
    fitOpts.measurements = generateMeasurements(line, t0, engine, true, true);
    fitOpts.startParameters = startParameters(line, fitOpts.measurements);
    fillPars(fitOpts.startParameters, recoY0, recoX0, recoTheta, recoPhi);
    //

    auto result = fitter.fit(std::move(fitOpts));
    if (!result.converged) {
      continue;
    }
    fillPars(result.parameters, recoY0, recoX0, recoTheta, recoPhi);
    covY0 = std::sqrt(result.covariance(toUnderlying(FitParIndex::y0),
                                        toUnderlying(FitParIndex::y0)));
    covX0 = std::sqrt(result.covariance(toUnderlying(FitParIndex::x0),
                                        toUnderlying(FitParIndex::x0)));
    covTheta = std::sqrt(result.covariance(toUnderlying(FitParIndex::theta),
                                           toUnderlying(FitParIndex::theta)));
    covPhi = std::sqrt(result.covariance(toUnderlying(FitParIndex::phi),
                                         toUnderlying(FitParIndex::phi)));

    chi2 = result.chi2;
    nDoF = result.nDoF;
    nIter = result.nIter;

    outTree->Fill();
  }
  outFile->WriteObject(outTree.get(), outTree->GetName());
  outTree.reset();
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Acts::Test

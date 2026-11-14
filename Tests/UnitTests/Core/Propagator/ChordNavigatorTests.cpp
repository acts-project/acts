// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/ChordNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;
using Acts::Experimental::ChordNavigator;

namespace ActsTests {

namespace {

const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
const MagneticFieldContext mctx{};

std::shared_ptr<const TrackingGeometry> gen3Geometry() {
  // The builder owns the detector elements of the surfaces
  static CylindricalTrackingGeometry builder(gctx, true);
  static const std::shared_ptr<const TrackingGeometry> geometry = [] {
    auto logger = getDefaultLogger("Geometry", Logging::ERROR);
    return std::shared_ptr<const TrackingGeometry>(builder(*logger));
  }();
  return geometry;
}

std::shared_ptr<const TrackingGeometry> gen3EndcapGeometry() {
  static CylindricalTrackingGeometry builder = [] {
    CylindricalTrackingGeometry b(gctx, true);
    b.gen3Endcaps = true;
    b.gen3CylinderPortalPolicy = true;
    return b;
  }();
  static const std::shared_ptr<const TrackingGeometry> geometry = [] {
    auto logger = getDefaultLogger("Geometry", Logging::ERROR);
    return std::shared_ptr<const TrackingGeometry>(builder(*logger));
  }();
  return geometry;
}

BoundTrackParameters startParameters(double phi, double theta, double qOverP) {
  return BoundTrackParameters::createCurvilinear(Vector4::Zero(), phi, theta,
                                                 qOverP, std::nullopt,
                                                 ParticleHypothesis::pion());
}

/// Sensitive surfaces reached by a propagation, with their positions
template <typename propagator_t>
std::vector<SurfaceHit> collectSensitive(const propagator_t& propagator,
                                         const BoundTrackParameters& start,
                                         double pathLimit,
                                         double maxStepSize = 0.) {
  using Actors = ActorList<SurfaceCollector<>>;
  typename propagator_t::template Options<Actors> options(gctx, mctx);
  options.pathLimit = pathLimit;
  options.loopProtection = false;
  options.maxSteps = 100000;
  if (maxStepSize > 0.) {
    options.stepping.maxStepSize = maxStepSize;
  }
  auto& collector = options.actorList.template get<SurfaceCollector<>>();
  collector.selector.selectSensitive = true;
  collector.selector.selectMaterial = false;
  auto result = propagator.propagate(start, options);
  BOOST_REQUIRE(result.ok());
  return result.value()
      .template get<SurfaceCollector<>::result_type>()
      .collected;
}

ChordNavigator makeChordNavigator(
    bool forwardHint,
    ChordNavigator::CandidateSearch search =
        ChordNavigator::CandidateSearch::Policy,
    bool envelopePreCheck = false,
    std::shared_ptr<const TrackingGeometry> geometry = nullptr) {
  ChordNavigator::Config cfg;
  cfg.trackingGeometry = geometry != nullptr ? geometry : gen3Geometry();
  cfg.candidateSearch = search;
  cfg.forwardHint = forwardHint;
  cfg.envelopePreCheck = envelopePreCheck;
  return ChordNavigator(cfg, getDefaultLogger("ChordNav", Logging::ERROR));
}

/// Whether @p hits holds a hit on the surface of @p hit at its position
bool contains(const std::vector<SurfaceHit>& hits, const SurfaceHit& hit) {
  return std::ranges::any_of(hits, [&hit](const SurfaceHit& candidate) {
    return candidate.surface == hit.surface &&
           (candidate.position - hit.position).norm() < 1_um;
  });
}

const std::vector<double> kPhis = {0.3, 1.4, 2.1, -2.6};
const std::vector<double> kThetas = {90_degree, 70_degree, 110_degree};

}  // namespace

BOOST_AUTO_TEST_SUITE(PropagatorSuite)

BOOST_AUTO_TEST_CASE(ChordNavigatorConstruction) {
  BOOST_CHECK_THROW(ChordNavigator(ChordNavigator::Config{}),
                    std::invalid_argument);

  CylindricalTrackingGeometry gen1Builder(gctx, false);
  ChordNavigator::Config gen1Cfg;
  gen1Cfg.trackingGeometry = gen1Builder();
  BOOST_CHECK_THROW(ChordNavigator{gen1Cfg}, std::invalid_argument);

  ChordNavigator::Config badCfg;
  badCfg.trackingGeometry = gen3Geometry();
  badCfg.minIncidence = 0.;
  BOOST_CHECK_THROW(ChordNavigator{badCfg}, std::invalid_argument);

  BOOST_CHECK_NO_THROW(makeChordNavigator(true));
}

// Without a field the chord is the trajectory: both navigators must reach the
// same surfaces at the same positions. The straight-line stepper has no step
// control. With the forward hint its steps end on the surfaces; without it,
// the step is bounded, because the crossings of the step that hits the path
// limit are not reported.
BOOST_AUTO_TEST_CASE(ChordNavigatorStraightLine) {
  const Navigator navigator({gen3Geometry()},
                            getDefaultLogger("Nav", Logging::ERROR));
  const Propagator<StraightLineStepper, Navigator> reference(
      StraightLineStepper{}, navigator);

  for (bool forwardHint : {true, false}) {
    const Propagator<StraightLineStepper, ChordNavigator> chord(
        StraightLineStepper{}, makeChordNavigator(forwardHint));
    const double maxStepSize = forwardHint ? 0. : 10_mm;
    for (double phi : kPhis) {
      for (double theta : kThetas) {
        const auto start = startParameters(phi, theta, 1. / 1_GeV);
        const auto expected = collectSensitive(reference, start, 1_m);
        const auto found = collectSensitive(chord, start, 1_m, maxStepSize);
        BOOST_TEST_CONTEXT("hint=" << forwardHint << " phi=" << phi
                                   << " theta=" << theta) {
          BOOST_CHECK(!expected.empty());
          BOOST_REQUIRE_EQUAL(found.size(), expected.size());
          for (std::size_t i = 0; i < found.size(); ++i) {
            BOOST_CHECK_EQUAL(found[i].surface, expected[i].surface);
            BOOST_CHECK_SMALL((found[i].position - expected[i].position).norm(),
                              1_um);
          }
        }
      }
    }
  }
}

// In a field, every crossing the straight-line prediction finds is also found
// on the chords, and every reached surface is hit inside its bounds. Crossings
// that coincide on a shared module edge can come in either order, so only
// containment is checked.
BOOST_AUTO_TEST_CASE(ChordNavigatorConstantField) {
  auto field = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));
  const Navigator navigator({gen3Geometry()},
                            getDefaultLogger("Nav", Logging::ERROR));
  const Propagator<EigenStepper<>, Navigator> reference(EigenStepper<>(field),
                                                        navigator);

  for (bool forwardHint : {true, false}) {
    const Propagator<EigenStepper<>, ChordNavigator> chord(
        EigenStepper<>(field), makeChordNavigator(forwardHint));
    for (double phi : kPhis) {
      for (double theta : kThetas) {
        for (double q : {-1., 1.}) {
          const auto start = startParameters(phi, theta, q / 1_GeV);
          const auto expected = collectSensitive(reference, start, 1_m);
          const auto found = collectSensitive(chord, start, 1_m);
          BOOST_TEST_CONTEXT("hint=" << forwardHint << " phi=" << phi
                                     << " theta=" << theta << " q=" << q) {
            BOOST_CHECK_GE(found.size(), expected.size());
            for (const SurfaceHit& hit : expected) {
              BOOST_CHECK_MESSAGE(contains(found, hit),
                                  "missing " << hit.surface->geometryId());
            }
            for (const SurfaceHit& hit : found) {
              BOOST_CHECK(hit.surface->isOnSurface(gctx, hit.position,
                                                   hit.direction,
                                                   BoundaryTolerance::None()));
            }
          }
        }
      }
    }
  }
}

// A 40 MeV track loops inside the barrel at 2 T. The straight-line prediction
// stops after the first half turn; the chords find the returning arms. The
// policy lookup must find the same crossings as testing every surface.
BOOST_AUTO_TEST_CASE(ChordNavigatorLooper) {
  auto field = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));

  const double pT = 40_MeV;
  const double theta = 87_degree;
  const double p = pT / std::sin(theta);
  // transverse radius of the circle and 3D length of one turn
  const double radius = pT / (0.299792458 * 2.) * 1_m / 1_GeV;
  const double turn = 2. * std::numbers::pi * radius / std::sin(theta);

  for (bool forwardHint : {true, false}) {
    const Propagator<EigenStepper<>, ChordNavigator> chord(
        EigenStepper<>(field), makeChordNavigator(forwardHint));
    const Propagator<EigenStepper<>, ChordNavigator> chordAll(
        EigenStepper<>(field),
        makeChordNavigator(forwardHint,
                           ChordNavigator::CandidateSearch::AllVolumeSurfaces));
    for (double phi : kPhis) {
      for (double q : {-1., 1.}) {
        const auto start = startParameters(phi, theta, q / p);
        const auto found = collectSensitive(chord, start, 2. * turn);
        const auto all = collectSensitive(chordAll, start, 2. * turn);
        BOOST_TEST_CONTEXT("hint=" << forwardHint << " phi=" << phi
                                   << " q=" << q) {
          // crossings beyond the first half turn
          std::size_t later = 0;
          for (const SurfaceHit& hit : found) {
            later += hit.position.z() / std::cos(theta) > 0.5 * turn ? 1 : 0;
            BOOST_CHECK(hit.surface->isOnSurface(
                gctx, hit.position, hit.direction, BoundaryTolerance::None()));
          }
          BOOST_CHECK_GT(later, 0u);
          BOOST_CHECK_EQUAL(found.size(), all.size());
          for (const SurfaceHit& hit : all) {
            BOOST_CHECK(contains(found, hit));
          }
        }
      }
    }
  }
}

// The envelope pre-check only skips chords that cannot cross a surface of the
// volume: loopers and straight tracks reach the same surfaces at the same
// positions with and without it.
BOOST_AUTO_TEST_CASE(ChordNavigatorEnvelopePreCheck) {
  ChordNavigator::Config badCfg;
  badCfg.trackingGeometry = gen3Geometry();
  badCfg.envelopePreCheck = true;
  badCfg.envelopeTolerance = -1.;
  BOOST_CHECK_THROW(ChordNavigator{badCfg}, std::invalid_argument);

  auto field = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));
  const ChordNavigator::CandidateSearch kPolicy =
      ChordNavigator::CandidateSearch::Policy;

  for (bool forwardHint : {true, false}) {
    const Propagator<EigenStepper<>, ChordNavigator> reference(
        EigenStepper<>(field), makeChordNavigator(forwardHint, kPolicy, false));
    const Propagator<EigenStepper<>, ChordNavigator> envelope(
        EigenStepper<>(field), makeChordNavigator(forwardHint, kPolicy, true));
    for (double pT : {25_MeV, 40_MeV, 1_GeV}) {
      for (double theta : {87_degree, 70_degree}) {
        for (double phi : kPhis) {
          for (double q : {-1., 1.}) {
            const auto start =
                startParameters(phi, theta, q * std::sin(theta) / pT);
            const auto expected = collectSensitive(reference, start, 1_m);
            const auto found = collectSensitive(envelope, start, 1_m);
            BOOST_TEST_CONTEXT("hint=" << forwardHint << " pT=" << pT
                                       << " theta=" << theta << " phi=" << phi
                                       << " q=" << q) {
              BOOST_CHECK(!expected.empty());
              BOOST_REQUIRE_EQUAL(found.size(), expected.size());
              for (std::size_t i = 0; i < found.size(); ++i) {
                BOOST_CHECK_EQUAL(found[i].surface, expected[i].surface);
                BOOST_CHECK_SMALL(
                    (found[i].position - expected[i].position).norm(), 1_um);
              }
            }
          }
        }
      }
    }
  }
}

// Forward loopers in the geometry with endcaps cross barrel layers and discs.
// The policy lookup, the reference search over all surfaces, and the envelope
// pre-check (with one policy query) find the same disc and barrel crossings.
BOOST_AUTO_TEST_CASE(ChordNavigatorEndcaps) {
  auto field = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));
  const auto geometry = gen3EndcapGeometry();
  const ChordNavigator::CandidateSearch kPolicy =
      ChordNavigator::CandidateSearch::Policy;
  const ChordNavigator::CandidateSearch kAll =
      ChordNavigator::CandidateSearch::AllVolumeSurfaces;

  for (bool forwardHint : {true, false}) {
    const Propagator<EigenStepper<>, ChordNavigator> all(
        EigenStepper<>(field),
        makeChordNavigator(forwardHint, kAll, false, geometry));
    const Propagator<EigenStepper<>, ChordNavigator> policy(
        EigenStepper<>(field),
        makeChordNavigator(forwardHint, kPolicy, false, geometry));
    ChordNavigator::Config fastCfg;
    fastCfg.trackingGeometry = geometry;
    fastCfg.forwardHint = forwardHint;
    fastCfg.queryFromChordEnd = false;
    fastCfg.envelopePreCheck = true;
    const Propagator<EigenStepper<>, ChordNavigator> fast(
        EigenStepper<>(field),
        ChordNavigator(fastCfg, getDefaultLogger("ChordNav", Logging::ERROR)));
    for (double pT : {25_MeV, 40_MeV}) {
      for (double cotTheta : {1.2, 3.}) {
        const double theta = std::atan2(1., cotTheta);
        for (double phi : kPhis) {
          for (double q : {-1., 1.}) {
            const auto start =
                startParameters(phi, theta, q * std::sin(theta) / pT);
            const auto expected = collectSensitive(all, start, 1_m);
            std::size_t discHits = 0;
            for (const SurfaceHit& hit : expected) {
              discHits += std::abs(hit.position.z()) > 550_mm ? 1 : 0;
            }
            BOOST_TEST_CONTEXT("hint=" << forwardHint << " pT=" << pT
                                       << " cot=" << cotTheta << " phi=" << phi
                                       << " q=" << q) {
              BOOST_CHECK_GT(discHits, 0u);
              for (const auto* propagator : {&policy, &fast}) {
                const auto found = collectSensitive(*propagator, start, 1_m);
                BOOST_REQUIRE_EQUAL(found.size(), expected.size());
                for (std::size_t i = 0; i < found.size(); ++i) {
                  BOOST_CHECK_EQUAL(found[i].surface, expected[i].surface);
                  BOOST_CHECK_SMALL(
                      (found[i].position - expected[i].position).norm(), 1_um);
                }
              }
            }
          }
        }
      }
    }
  }
}

/// Moves the state within the first sensitive surface reached, like a Kalman
/// update with a large residual, and keeps direction and momentum
struct InPlaneJump {
  struct this_result {
    bool jumped = false;
    Vector3 position = Vector3::Zero();
    Vector3 direction = Vector3::Zero();
    double qOverP = 0.;
    double pathLength = 0.;
    const Surface* surface = nullptr;
  };
  using result_type = this_result;

  /// Shift along the local x and y axes of the surface
  double alongX = 0.;
  double alongY = 0.;

  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  Result<void> act(propagator_state_t& state, const stepper_t& stepper,
                   const navigator_t& navigator, result_type& result,
                   const Logger& /*logger*/) const {
    if (result.jumped || state.stage == PropagatorStage::postPropagation) {
      return {};
    }
    const Surface* surface = navigator.currentSurface(state.navigation);
    if (surface == nullptr || !surface->isSensitive()) {
      return {};
    }
    const auto& rotation =
        surface->localToGlobalTransform(state.geoContext).linear();
    const Vector3 position = stepper.position(state.stepping) +
                             alongX * rotation.col(0) +
                             alongY * rotation.col(1);
    const Vector3 direction = stepper.direction(state.stepping);
    stepper.update(state.stepping, position, direction,
                   stepper.qOverP(state.stepping),
                   stepper.time(state.stepping));
    result = {true,
              position,
              direction,
              stepper.qOverP(state.stepping),
              stepper.pathLength(state.stepping),
              surface};
    return {};
  }
};

// An actor moves the state within a reached surface (as a Kalman update
// does). The chord search must not span that jump: after it, the navigator
// finds exactly the hits of a propagation started from the moved state.
BOOST_AUTO_TEST_CASE(ChordNavigatorActorJump) {
  auto field = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));
  const double pT = 40_MeV;
  const double theta = 87_degree;
  const double p = pT / std::sin(theta);
  const double radius = pT / (0.299792458 * 2.) * 1_m / 1_GeV;
  const double turn = 2. * std::numbers::pi * radius / std::sin(theta);
  const double pathLimit = 1.5 * turn;

  using Actors = ActorList<SurfaceCollector<>, InPlaneJump>;
  std::size_t nJumps = 0;
  for (bool forwardHint : {true, false}) {
    const Propagator<EigenStepper<>, ChordNavigator> chord(
        EigenStepper<>(field), makeChordNavigator(forwardHint));
    for (double phi : kPhis) {
      for (double q : {-1., 1.}) {
        for (const auto& [alongX, alongY] :
             std::vector<std::pair<double, double>>{
                 {-12_mm, 0.}, {-6_mm, 15_mm}, {6_mm, -15_mm}, {12_mm, 0.}}) {
          BOOST_TEST_CONTEXT("hint=" << forwardHint << " phi=" << phi
                                     << " q=" << q << " jump=(" << alongX << ","
                                     << alongY << ")") {
            Propagator<EigenStepper<>, ChordNavigator>::Options<Actors> options(
                gctx, mctx);
            options.pathLimit = pathLimit;
            options.loopProtection = false;
            options.maxSteps = 100000;
            auto& collector =
                options.actorList.template get<SurfaceCollector<>>();
            collector.selector.selectSensitive = true;
            collector.selector.selectMaterial = false;
            auto& jump = options.actorList.template get<InPlaneJump>();
            jump.alongX = alongX;
            jump.alongY = alongY;

            const auto start = startParameters(phi, theta, q / p);
            auto result = chord.propagate(start, options);
            BOOST_REQUIRE(result.ok());
            const auto& jumped = result.value().get<InPlaneJump::result_type>();
            BOOST_REQUIRE(jumped.jumped);
            ++nJumps;
            const auto& hits =
                result.value().get<SurfaceCollector<>::result_type>().collected;
            // the hits after the one where the state jumped
            const auto at = std::ranges::find_if(hits, [&](const auto& hit) {
              return hit.surface == jumped.surface;
            });
            BOOST_REQUIRE(at != hits.end());
            const std::vector<SurfaceHit> after(std::next(at), hits.end());

            // reference: a fresh propagation from the moved state
            const auto restart = BoundTrackParameters::createCurvilinear(
                Vector4(jumped.position.x(), jumped.position.y(),
                        jumped.position.z(), 0.),
                jumped.direction, jumped.qOverP, std::nullopt,
                ParticleHypothesis::pion());
            const auto expected =
                collectSensitive(chord, restart, pathLimit - jumped.pathLength);
            BOOST_CHECK_EQUAL(after.size(), expected.size());
            for (const SurfaceHit& hit : expected) {
              BOOST_CHECK(contains(after, hit));
            }
          }
        }
      }
    }
  }
  BOOST_CHECK_EQUAL(nJumps, 2u * kPhis.size() * 2u * 4u);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

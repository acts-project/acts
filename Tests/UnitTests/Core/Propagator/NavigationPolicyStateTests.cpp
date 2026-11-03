// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldError.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "ActsTests/CommonHelpers/TemporaryDirectory.hpp"

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <format>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;
using Acts::Blueprint;
using Acts::BlueprintOptions;
using Acts::StaticBlueprintNode;

namespace ActsTests {

GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
MagneticFieldContext mfContext = MagneticFieldContext();

void step(Vector3& pos, const Vector3& dir, const Surface& surface) {
  Intersection3D intersection =
      surface.intersect(gctx, pos, dir).closestForward();
  pos += intersection.pathLength() * dir;
}

void step(Vector3& pos, const Vector3& dir, const NavigationTarget& target) {
  step(pos, dir, target.surface());
}

// ---------------------------------------------------------------------------
// Custom Policy 1: AlternatingInvalidPolicy
//
// Returns isValid()=false every other call. Does NOT add navigation
// candidates -- composed with TryAllNavigationPolicy via the factory.
// ---------------------------------------------------------------------------

class AlternatingInvalidPolicy final : public INavigationPolicy {
 public:
  struct State {
    int counter = 0;
  };

  AlternatingInvalidPolicy(const GeometryContext& /*gctx*/,
                           const TrackingVolume& /*volume*/,
                           const Logger& /*logger*/) {}

  void initializeCandidates(const GeometryContext& /*gctx*/,
                            const NavigationArguments& /*args*/,
                            NavigationPolicyState& /*state*/,
                            AppendOnlyNavigationStream& /*stream*/,
                            const Logger& /*logger*/) const {
    // No-op: TryAllNavigationPolicy provides the candidates.
  }

  void connect(NavigationDelegate& delegate) const override {
    connectDefault<AlternatingInvalidPolicy>(delegate);
  }

  bool isValid(const GeometryContext& /*gctx*/,
               const NavigationArguments& /*args*/,
               NavigationPolicyState& state,
               const Logger& logger) const override {
    auto& s = state.as<State>();
    s.counter++;
    bool valid = (s.counter % 2 == 0);
    ACTS_VERBOSE("AlternatingInvalidPolicy isValid: counter="
                 << s.counter << " valid=" << valid);
    return valid;
  }

  void createState(const GeometryContext& /*gctx*/,
                   const NavigationArguments& /*args*/,
                   NavigationPolicyStateManager& stateManager,
                   const Logger& logger) const override {
    ACTS_VERBOSE("AlternatingInvalidPolicy createState");
    stateManager.pushState<State>();
  }

  void popState(NavigationPolicyStateManager& stateManager,
                const Logger& logger) const override {
    ACTS_VERBOSE("AlternatingInvalidPolicy popState");
    stateManager.popState();
  }
};

static_assert(NavigationPolicyConcept<AlternatingInvalidPolicy>);

// ---------------------------------------------------------------------------
// Custom Policy 2: ConeValidityPolicy
//
// Stores the direction at createState time. isValid returns false when
// the current direction deviates from the stored direction by more than
// a configurable cone angle (dot-product check).
// ---------------------------------------------------------------------------

class ConeValidityPolicy final : public INavigationPolicy {
 public:
  struct State {
    Vector3 initialDirection = Vector3::Zero();
  };

  ConeValidityPolicy(const GeometryContext& /*gctx*/,
                     const TrackingVolume& /*volume*/, const Logger& /*logger*/,
                     double coneAngle,
                     std::shared_ptr<std::atomic<int>> reresolutionCounter)
      : m_cosConeAngle(std::cos(coneAngle)),
        m_reresolutionCounter(std::move(reresolutionCounter)) {}

  void initializeCandidates(const GeometryContext& /*gctx*/,
                            const NavigationArguments& args,
                            NavigationPolicyState& state,
                            AppendOnlyNavigationStream& /*stream*/,
                            const Logger& logger) const {
    // Re-center the cone on the current direction after re-resolution.
    auto& s = state.as<State>();
    s.initialDirection = args.direction.normalized();
    ACTS_VERBOSE("ConeValidityPolicy initializeCandidates: re-centered to "
                 << s.initialDirection.transpose());
  }

  void connect(NavigationDelegate& delegate) const override {
    connectDefault<ConeValidityPolicy>(delegate);
  }

  bool isValid(const GeometryContext& /*gctx*/, const NavigationArguments& args,
               NavigationPolicyState& state,
               const Logger& logger) const override {
    auto& s = state.as<State>();
    double dot = s.initialDirection.dot(args.direction.normalized());
    bool valid = (dot >= m_cosConeAngle);
    if (!valid) {
      ++(*m_reresolutionCounter);
      ACTS_VERBOSE("ConeValidityPolicy isValid: INVALID dot="
                   << dot << " threshold=" << m_cosConeAngle);
    } else {
      ACTS_VERBOSE("ConeValidityPolicy isValid: VALID dot="
                   << dot << " threshold=" << m_cosConeAngle);
    }
    return valid;
  }

  void createState(const GeometryContext& /*gctx*/,
                   const NavigationArguments& args,
                   NavigationPolicyStateManager& stateManager,
                   const Logger& logger) const override {
    ACTS_VERBOSE("ConeValidityPolicy createState: direction="
                 << args.direction.transpose());
    auto& s = stateManager.pushState<State>();
    s.initialDirection = args.direction.normalized();
  }

  void popState(NavigationPolicyStateManager& stateManager,
                const Logger& logger) const override {
    ACTS_VERBOSE("ConeValidityPolicy popState");
    stateManager.popState();
  }

 private:
  double m_cosConeAngle;
  std::shared_ptr<std::atomic<int>> m_reresolutionCounter;
};

static_assert(NavigationPolicyConcept<ConeValidityPolicy>);

// ---------------------------------------------------------------------------
// StepCollector actor (records trajectory positions for visualization)
// ---------------------------------------------------------------------------

struct StepCollector {
  struct this_result {
    std::vector<Vector3> position;
    std::vector<const Surface*> surfaces;
  };

  using result_type = this_result;

  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  Result<void> act(propagator_state_t& state, const stepper_t& stepper,
                   const navigator_t& /*navigator*/, result_type& result,
                   const Logger& /*logger*/) const {
    result.position.push_back(stepper.position(state.stepping));
    if (state.navigation.currentSurface != nullptr) {
      result.surfaces.push_back(state.navigation.currentSurface);
    }
    return Result<void>::success();
  }

  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  bool checkAbort(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                  const navigator_t& /*navigator*/, result_type& /*result*/,
                  const Logger& /*logger*/) const {
    return false;
  }
};

// ---------------------------------------------------------------------------
// Geometry builder: outer cuboid with three inner cuboids
// ---------------------------------------------------------------------------

std::unique_ptr<TrackingGeometry> buildBoxGeometry(const BlueprintOptions& opts,
                                                   const Logger& logger) {
  Blueprint::Config cfg;
  cfg.envelope = ExtentEnvelope{{
      .x = {20_mm, 20_mm},
      .y = {20_mm, 20_mm},
      .z = {20_mm, 20_mm},
  }};

  Blueprint root{cfg};

  // Outer parent volume: 2m x 2m x 2m at origin
  auto outerBounds = std::make_shared<CuboidVolumeBounds>(1_m, 1_m, 1_m);
  auto outerVol = std::make_unique<TrackingVolume>(Transform3::Identity(),
                                                   outerBounds, "Outer");

  auto outerNode = std::make_shared<StaticBlueprintNode>(std::move(outerVol));

  // Box1: left side along -X
  auto box1Bounds =
      std::make_shared<CuboidVolumeBounds>(300_mm, 300_mm, 300_mm);
  auto box1Vol = std::make_unique<TrackingVolume>(
      Transform3(Translation3(-500_mm, 0, 0)), box1Bounds, "BoxLeft");

  auto box1Node = std::make_shared<StaticBlueprintNode>(std::move(box1Vol));

  // Box2: right side along +X, centered vertically
  auto box2Bounds =
      std::make_shared<CuboidVolumeBounds>(300_mm, 300_mm, 300_mm);
  auto box2Vol = std::make_unique<TrackingVolume>(
      Transform3(Translation3(500_mm, 0, 0)), box2Bounds, "BoxRight");

  auto box2Node = std::make_shared<StaticBlueprintNode>(std::move(box2Vol));

  // Box3: above Box2, shifted in +Y
  auto box3Bounds =
      std::make_shared<CuboidVolumeBounds>(200_mm, 200_mm, 200_mm);
  auto box3Vol = std::make_unique<TrackingVolume>(
      Transform3(Translation3(500_mm, 600_mm, 0)), box3Bounds, "BoxTop");

  auto box3Node = std::make_shared<StaticBlueprintNode>(std::move(box3Vol));

  outerNode->addChild(box1Node);
  outerNode->addChild(box2Node);
  outerNode->addChild(box3Node);

  root.addChild(outerNode);

  return root.construct(opts, gctx, logger);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NavigationPolicyStateSuite)

// ===== Test 1: AlternatingInvalidPolicy with straight-line stepping =========

BOOST_AUTO_TEST_CASE(AlternatingInvalidPolicy_StraightLine) {
  TemporaryDirectory tmp{};
  auto logger = getDefaultLogger("AlternatingTest", Logging::VERBOSE);

  // Build geometry with AlternatingInvalidPolicy on all volumes
  BlueprintOptions opts;
  opts.defaultNavigationPolicyFactory = NavigationPolicyFactory{}
                                            .add<TryAllNavigationPolicy>()
                                            .add<AlternatingInvalidPolicy>()
                                            .asUniquePtr();

  auto trackingGeometry = buildBoxGeometry(opts, *logger);
  BOOST_REQUIRE(trackingGeometry);
  BOOST_CHECK(trackingGeometry->geometryVersion() ==
              TrackingGeometry::GeometryVersion::Gen3);

  // Write OBJ output
  ObjVisualization3D vis;
  trackingGeometry->visualize(vis, gctx);
  vis.write(tmp.path() / "alternating_invalid_policy_geometry.obj");

  // Set up Navigator
  auto sharedGeometry =
      std::shared_ptr<const TrackingGeometry>(std::move(trackingGeometry));

  Navigator::Config navCfg;
  navCfg.trackingGeometry = sharedGeometry;
  navCfg.resolveSensitive = true;
  navCfg.resolveMaterial = true;
  navCfg.resolvePassive = false;
  Navigator navigator(navCfg, logger->clone("Navigator"));

  // Start at left side of outer volume, propagate along +X
  Navigator::Options options(gctx);
  Navigator::State state = navigator.makeState(options);

  Vector3 position{-900_mm, 0, 0};
  Vector3 direction = Vector3::UnitX();

  Result<void> result = navigator.initialize(
      state, {.position = position,
              .direction = direction,
              .propagationDirection = Direction::Forward()});
  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE(state.currentVolume != nullptr);

  // Step through geometry, track visited volumes
  std::vector<std::string> visitedVolumes;
  visitedVolumes.push_back(state.currentVolume->volumeName());

  int maxSteps = 200;
  for (int i = 0; i < maxSteps; ++i) {
    NavigationTarget target = navigator.nextTarget(state, position, direction);

    if (target.isNone()) {
      break;
    }

    step(position, direction, target);
    navigator.handleSurfaceReached(state, position, direction,
                                   target.surface());

    if (state.currentVolume != nullptr &&
        (visitedVolumes.empty() ||
         visitedVolumes.back() != state.currentVolume->volumeName())) {
      visitedVolumes.push_back(state.currentVolume->volumeName());
    }

    if (state.navigationBreak) {
      break;
    }
  }

  // Navigation should have completed (left the world)
  BOOST_CHECK(state.navigationBreak || state.currentVolume == nullptr);

  // Should have visited at least: Outer -> BoxLeft -> Outer -> BoxRight ->
  // Outer
  BOOST_CHECK_GE(visitedVolumes.size(), 3u);

  // Print visited volumes for diagnostics
  for (const auto& name : visitedVolumes) {
    BOOST_TEST_MESSAGE("Visited: " << name);
  }
}

// ===== Test 2: ConeValidityPolicy with magnetic field ======================

BOOST_AUTO_TEST_CASE(ConeValidityPolicy_MagneticField) {
  TemporaryDirectory tmp{};
  auto logger = getDefaultLogger("ConeTest", Logging::VERBOSE);

  const double coneAngle = 10.0 * std::numbers::pi / 180.0;

  using EigenStepperType = EigenStepper<>;
  using EigenPropagatorType = Propagator<EigenStepperType, Navigator>;

  // Helper lambda to run a single sub-case
  auto runSubCase = [&](double bFieldZ, const std::string& label) -> int {
    auto counter = std::make_shared<std::atomic<int>>(0);

    BlueprintOptions opts;
    opts.defaultNavigationPolicyFactory =
        NavigationPolicyFactory{}
            .add<TryAllNavigationPolicy>()
            .add<ConeValidityPolicy>(coneAngle, counter)
            .asUniquePtr();

    auto trackingGeometry = buildBoxGeometry(opts, *logger);
    BOOST_REQUIRE(trackingGeometry);

    ObjVisualization3D vis;
    trackingGeometry->visualize(vis, gctx);
    vis.write(tmp.path() / "cone_validity_geometry.obj");

    auto sharedGeometry =
        std::shared_ptr<const TrackingGeometry>(std::move(trackingGeometry));

    auto bField = std::make_shared<ConstantBField>(Vector3{0, 0, bFieldZ});
    EigenStepperType stepper(bField);

    Navigator::Config navCfg;
    navCfg.trackingGeometry = sharedGeometry;
    navCfg.resolveSensitive = true;
    navCfg.resolveMaterial = true;
    navCfg.resolvePassive = false;
    Navigator navigator(navCfg, logger->clone("Navigator"));

    EigenPropagatorType propagator(std::move(stepper), std::move(navigator),
                                   logger->clone("Propagator"));

    // 1 GeV pion starting at left side, direction +X
    Vector3 startPos{-900_mm, 0, 0};
    Vector3 startDir = Vector3::UnitX();

    BoundTrackParameters start = BoundTrackParameters::createCurvilinear(
        Vector4(startPos.x(), startPos.y(), startPos.z(), 0), startDir,
        1.0 / 1_GeV, std::nullopt, ParticleHypothesis::pion());

    using MyActorList = ActorList<StepCollector, EndOfWorldReached>;
    EigenPropagatorType::Options<MyActorList> propOpts(gctx, mfContext);
    propOpts.pathLimit = 5_m;
    propOpts.stepping.maxStepSize = 100_mm;

    auto result = propagator.propagate(start, propOpts);
    BOOST_REQUIRE(result.ok());

    const auto& info =
        result.value().template get<StepCollector::this_result>();

    const auto& positions = info.position;
    if (!positions.empty()) {
      ObjVisualization3D trajVis;
      for (std::size_t i = 0; i + 1 < positions.size(); ++i) {
        trajVis.line(positions[i], positions[i + 1]);
      }
      trajVis.write(tmp.path() / std::format("trajectory_{}.obj", label));
    }

    ObjVisualization3D srfVis;
    if (!info.surfaces.empty()) {
      for (const auto* srf : info.surfaces) {
        std::cout << srf->toString(gctx) << std::endl;
        if (srf->type() != Surface::Plane ||
            srf->bounds().type() != SurfaceBounds::eRectangle) {
          continue;
        }
        srf->visualize(srfVis, gctx);
      }
      srfVis.write(tmp.path() / std::format("surfaces_{}.obj", label));
    }

    return counter->load();
  };

  // Sub-case a: strong field
  int strongReresolutions = runSubCase(-1.6_T, "strong");
  BOOST_TEST_MESSAGE("Strong field re-resolutions: " << strongReresolutions);
  BOOST_CHECK_GT(strongReresolutions, 0);

  // Sub-case b: weak field
  int weakReresolutions = runSubCase(0.001_T, "weak");
  BOOST_TEST_MESSAGE("Weak field re-resolutions: " << weakReresolutions);

  // Strong field must cause more re-resolutions than weak field
  BOOST_CHECK_GT(strongReresolutions, weakReresolutions);
}

// Records all three policy entry points, including entry through portals and
// forwarding through MultiNavigationPolicy. The synthetic construction-time
// probes are discarded before propagation.
struct CurvatureRecord {
  std::string entry;
  std::string volume;
  NavigationArguments args;
  double qOverP;
};

struct CurvatureRecords {
  std::vector<CurvatureRecord> calls;
  double qOverP = 1. / 1_GeV;
};

// Hide the curvature overloads to model an existing custom navigator.
class LegacyNavigator : public DirectNavigator {
 public:
  struct Calls {
    unsigned int next = 0;
    unsigned int valid = 0;
    unsigned int reached = 0;
  };
  std::shared_ptr<Calls> calls = std::make_shared<Calls>();

  NavigationTarget nextTarget(State& state, const Vector3& position,
                              const Vector3& direction) const {
    ++calls->next;
    return DirectNavigator::nextTarget(state, position, direction);
  }
  bool checkTargetValid(const State& state, const Vector3& position,
                        const Vector3& direction) const {
    ++calls->valid;
    return DirectNavigator::checkTargetValid(state, position, direction);
  }
  void handleSurfaceReached(State& state, const Vector3& position,
                            const Vector3& direction,
                            const Surface& surface) const {
    ++calls->reached;
    DirectNavigator::handleSurfaceReached(state, position, direction, surface);
  }
};

static_assert(NavigatorConcept<LegacyNavigator>);
static_assert(NavigatorConcept<Navigator>);

BOOST_AUTO_TEST_CASE(LegacyNavigatorPropagation) {
  LegacyNavigator navigator;
  const auto calls = navigator.calls;
  Propagator propagator(StraightLineStepper{}, std::move(navigator));
  const auto surface =
      CurvilinearSurface(Vector3(0., 0., 10_mm), Vector3::UnitZ())
          .planeSurface();
  decltype(propagator)::Options<> options(gctx, mfContext);
  options.navigation.surfaceSequence = {surface.get()};
  options.stepping.maxStepSize = 1_mm;
  options.pathLimit = 20_mm;
  const auto start = BoundTrackParameters::createCurvilinear(
      Vector4::Zero(), Vector3::UnitZ(), 1. / 1_GeV, std::nullopt,
      ParticleHypothesis::pion());
  const auto result = propagator.propagate(start, options);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_GT(calls->next, 0u);
  BOOST_CHECK_GT(calls->valid, 0u);
  BOOST_CHECK_EQUAL(calls->reached, 1u);
}

class CurvatureRecordingPolicy final : public INavigationPolicy {
 public:
  struct State {};

  CurvatureRecordingPolicy(const GeometryContext& /*gctx*/,
                           const TrackingVolume& volume,
                           const Logger& /*logger*/,
                           std::shared_ptr<CurvatureRecords> records)
      : m_volume(volume.volumeName()), m_records(std::move(records)) {}

  void connect(NavigationDelegate& delegate) const override {
    connectDefault<CurvatureRecordingPolicy>(delegate);
  }

  void createState(const GeometryContext& /*gctx*/,
                   const NavigationArguments& args,
                   NavigationPolicyStateManager& manager,
                   const Logger& /*logger*/) const override {
    record("create", args);
    manager.pushState<State>();
  }

  bool isValid(const GeometryContext& /*gctx*/, const NavigationArguments& args,
               NavigationPolicyState& /*state*/,
               const Logger& /*logger*/) const override {
    record("valid", args);
    return true;
  }

  void initializeCandidates(const GeometryContext& /*gctx*/,
                            const NavigationArguments& args,
                            NavigationPolicyState& /*state*/,
                            AppendOnlyNavigationStream& /*stream*/,
                            const Logger& /*logger*/) const {
    record("candidates", args);
  }

 private:
  void record(const std::string& entry, const NavigationArguments& args) const {
    m_records->calls.push_back({entry, m_volume, args, m_records->qOverP});
  }
  std::string m_volume;
  std::shared_ptr<CurvatureRecords> m_records;
};

// Change both the momentum and the direction before the first target query,
// then change the momentum again after the first step. Policy callbacks must
// see the updated curvature, not a derivative cached by the previous step.
struct ChangeCurvatureActor {
  struct result_type {
    bool changedAfterStep = false;
  };
  std::shared_ptr<CurvatureRecords> records;

  template <typename PropagatorState, typename Stepper, typename Nav>
  Result<void> act(PropagatorState& state, const Stepper& stepper,
                   const Nav& /*navigator*/, result_type& result,
                   const Logger& /*logger*/) const {
    Vector3 direction = stepper.direction(state.stepping);
    if (state.stage == PropagatorStage::prePropagation) {
      direction = Vector3(1., 0.02, 0.01).normalized();
    } else if (state.stage == PropagatorStage::postStep &&
               !result.changedAfterStep) {
      result.changedAfterStep = true;
    } else {
      return Result<void>::success();
    }
    records->qOverP *= 0.5;
    stepper.update(state.stepping, stepper.position(state.stepping), direction,
                   records->qOverP, stepper.time(state.stepping));
    return Result<void>::success();
  }

  template <typename PropagatorState, typename Stepper, typename Nav>
  bool checkAbort(PropagatorState& /*state*/, const Stepper& /*stepper*/,
                  const Nav& /*navigator*/, result_type& /*result*/,
                  const Logger& /*logger*/) const {
    return false;
  }
};

BOOST_AUTO_TEST_CASE(CurvatureReachesEveryPolicyEntryPoint) {
  for (Direction propagationDirection :
       {Direction::Forward(), Direction::Backward()}) {
    for (double chargeSign : {-1., 1.}) {
      auto records = std::make_shared<CurvatureRecords>();
      records->qOverP = chargeSign / 1_GeV;
      BlueprintOptions opts;
      opts.defaultNavigationPolicyFactory =
          NavigationPolicyFactory{}
              .add<TryAllNavigationPolicy>()
              .add<CurvatureRecordingPolicy>(records)
              .asUniquePtr();
      std::shared_ptr<const Logger> logger =
          getDefaultLogger("CurvatureTest", Logging::WARNING);
      std::shared_ptr<const TrackingGeometry> geometry =
          buildBoxGeometry(opts, *logger);
      records->calls.clear();
      const Vector3 field(0., 0., 0.1_T);
      Navigator navigator({.trackingGeometry = geometry}, logger);
      Propagator propagator(
          EigenStepper<>(std::make_shared<ConstantBField>(field)),
          std::move(navigator), logger);
      const auto start = BoundTrackParameters::createCurvilinear(
          Vector4(
              propagationDirection == Direction::Forward() ? -900_mm : 900_mm,
              0., 0., 0.),
          Vector3::UnitX(), records->qOverP, std::nullopt,
          ParticleHypothesis::pion());
      decltype(propagator)::Options<
          ActorList<ChangeCurvatureActor, EndOfWorldReached>>
          options(gctx, mfContext);
      options.direction = propagationDirection;
      options.pathLimit = 3_m;
      options.stepping.maxStepSize = 50_mm;
      options.actorList.get<ChangeCurvatureActor>().records = records;
      auto result = propagator.propagate(start, options);
      BOOST_REQUIRE(result.ok());
      BOOST_REQUIRE(!records->calls.empty());
      BOOST_CHECK_SMALL((records->calls.front().args.curvature -
                         Vector3(0., -chargeSign * 0.1_T / 1_GeV, 0.))
                            .norm(),
                        1e-12);
      std::size_t creates = 0, candidates = 0, valid = 0;
      bool sawMomentumUpdate = false;
      bool sawNewVolume = false;
      for (const auto& call : records->calls) {
        const Vector3 physicalDirection =
            propagationDirection * call.args.direction;
        const Vector3 expected = call.qOverP * physicalDirection.cross(field);
        BOOST_CHECK_SMALL((call.args.curvature - expected).norm(), 1e-12);
        BOOST_CHECK_SMALL(call.args.curvature.dot(call.args.direction), 1e-12);
        if (call.entry == "create") {
          ++creates;
        }
        if (call.entry == "candidates") {
          ++candidates;
        }
        if (call.entry == "valid") {
          ++valid;
        }
        sawMomentumUpdate |= std::abs(call.qOverP) < 0.3 / 1_GeV;
        sawNewVolume |= call.entry == "create" &&
                        call.volume != records->calls.front().volume;
      }
      BOOST_CHECK_GT(creates, 1u);
      BOOST_CHECK_GT(candidates, 0u);
      BOOST_CHECK_GT(valid, 0u);
      BOOST_CHECK(sawMomentumUpdate);
      BOOST_CHECK(sawNewVolume);
    }
  }
}

class FailingCurvatureField final : public MagneticFieldProvider {
 public:
  bool fail = true;
  Cache makeCache(const MagneticFieldContext& /*mctx*/) const override {
    return Cache{};
  }
  Result<Vector3> getField(const Vector3& /*position*/,
                           Cache& /*cache*/) const override {
    if (fail) {
      return Result<Vector3>::failure(MagneticFieldError::OutOfBounds);
    }
    return Result<Vector3>::success(Vector3::Zero());
  }
};

// Model a custom stepper written before the curvature accessor existed.
class LegacyAtlasStepper : public AtlasStepper {
 public:
  using AtlasStepper::AtlasStepper;
  Result<Vector3> curvature(State& state) const = delete;
};
static_assert(StepperConcept<LegacyAtlasStepper>);

// An accessor must take precedence over the compatibility field calculation.
class CustomCurvatureStepper : public AtlasStepper {
 public:
  using AtlasStepper::AtlasStepper;
  Result<Vector3> curvature(State& /*state*/) const {
    return Result<Vector3>::success(Vector3(0., 0.25, 0.));
  }
};
static_assert(StepperConcept<CustomCurvatureStepper>);

BOOST_AUTO_TEST_CASE(StepperCurvatureCompatibility) {
  const auto start = BoundTrackParameters::createCurvilinear(
      Vector4::Zero(), Vector3::UnitX(), 1. / 1_GeV, std::nullopt,
      ParticleHypothesis::pion());
  const Vector3 field(0., 0., 1_T);
  Propagator legacy{
      LegacyAtlasStepper(std::make_shared<ConstantBField>(field))};
  decltype(legacy)::Options<> options(gctx, mfContext);
  auto state = legacy.makeState(options);
  BOOST_REQUIRE(legacy.initialize(state, start, nullptr).ok());
  BOOST_CHECK_SMALL((state.curvature - Vector3(0., -1_T / 1_GeV, 0.)).norm(),
                    1e-18);

  Propagator custom{
      CustomCurvatureStepper(std::make_shared<FailingCurvatureField>())};
  decltype(custom)::Options<> customOptions(gctx, mfContext);
  customOptions.loopProtection = false;
  auto customState = custom.makeState(customOptions);
  BOOST_REQUIRE(custom.initialize(customState, start, nullptr).ok());
  BOOST_CHECK_SMALL((customState.curvature - Vector3(0., 0.25, 0.)).norm(),
                    1e-18);
}

BOOST_AUTO_TEST_CASE(CurvatureZeroCasesAndFieldErrors) {
  const auto start = [](ParticleHypothesis particle =
                            ParticleHypothesis::pion()) {
    return BoundTrackParameters::createCurvilinear(
        Vector4::Zero(), Vector3::UnitX(),
        particle.qOverP(1_GeV, particle.absoluteCharge()), std::nullopt,
        particle);
  };
  // Exercise the propagator initialization separately: neutral EigenStepper
  // integration is not part of this test, only the supplied navigation data.
  const auto checkZero = [&](auto stepper, ParticleHypothesis particle) {
    Propagator propagator(std::move(stepper));
    typename decltype(propagator)::template Options<> options(gctx, mfContext);
    options.loopProtection = false;
    auto state = propagator.makeState(options);
    state.curvature = Vector3::Ones();  // Must not leak across initialization.
    BOOST_REQUIRE(propagator.initialize(state, start(particle), nullptr).ok());
    // Bound angle conversion can leave a roundoff-sized component transverse
    // to a nominally parallel field.
    BOOST_CHECK_SMALL(state.curvature.norm(), 1e-18);
  };
  checkZero(StraightLineStepper{}, ParticleHypothesis::pion());
  checkZero(EigenStepper<>(std::make_shared<ConstantBField>(Vector3::Zero())),
            ParticleHypothesis::pion());
  checkZero(
      EigenStepper<>(std::make_shared<ConstantBField>(Vector3(1_T, 0., 0.))),
      ParticleHypothesis::pion());
  checkZero(EigenStepper<>(std::make_shared<FailingCurvatureField>()),
            ParticleHypothesis::photon().withMomentumHypothesis(1_GeV));

  Propagator propagator(
      EigenStepper<>(std::make_shared<FailingCurvatureField>()));
  decltype(propagator)::Options<> options(gctx, mfContext);
  auto state = propagator.makeState(options);
  const auto result = propagator.initialize(state, start(), nullptr);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() ==
              make_error_code(MagneticFieldError::OutOfBounds));

  // A field lookup can also fail when refreshing an already initialized run.
  auto field = std::make_shared<FailingCurvatureField>();
  field->fail = false;
  Propagator refreshingPropagator{EigenStepper<>(field)};
  auto refreshingState = refreshingPropagator.makeState(options);
  BOOST_REQUIRE(
      refreshingPropagator.initialize(refreshingState, start(), nullptr).ok());
  field->fail = true;
  const auto refreshResult = refreshingPropagator.propagate(refreshingState);
  BOOST_REQUIRE(!refreshResult.ok());
  BOOST_CHECK(refreshResult.error() ==
              make_error_code(MagneticFieldError::OutOfBounds));
}

// Custom steppers with the original surface-status signature remain usable.
class LegacySurfaceStatusStepper : public AtlasStepper {
 public:
  using AtlasStepper::AtlasStepper;
  IntersectionStatus updateSurfaceStatus(
      State& state, const Surface& surface, std::uint8_t index,
      Direction direction, const BoundaryTolerance& bounds, double tolerance,
      ConstrainedStep::Type type, const Logger& logger) const {
    return AtlasStepper::updateSurfaceStatus(state, surface, index, direction,
                                             bounds, tolerance, type, logger);
  }
};
static_assert(StepperConcept<LegacySurfaceStatusStepper>);

// Observe the dispatched curvature while deliberately retaining a ray check.
class SurfaceStatusCurvatureRecorder : public AtlasStepper {
 public:
  using AtlasStepper::AtlasStepper;
  mutable Vector3 receivedCurvature = Vector3::Zero();
  IntersectionStatus updateSurfaceStatus(
      State& state, const Surface& surface, std::uint8_t index,
      Direction direction, const BoundaryTolerance& bounds, double tolerance,
      ConstrainedStep::Type type, const Logger& logger,
      const Vector3& curvature) const {
    receivedCurvature = curvature;
    return AtlasStepper::updateSurfaceStatus(state, surface, index, direction,
                                             bounds, tolerance, type, logger);
  }
};

BOOST_AUTO_TEST_CASE(CurvedTargetStatusCompatibilityAndErrors) {
  const auto start = BoundTrackParameters::createCurvilinear(
      Vector4::Zero(), Vector3::UnitZ(), 1. / 1_GeV, std::nullopt,
      ParticleHypothesis::pion());
  const auto target = Surface::makeShared<PlaneSurface>(
      Transform3(Translation3(0., 1., 20.)),
      std::make_shared<RectangleBounds>(0.1, 0.1));
  const auto check = [&](const auto& stepper, auto& state) {
    return detail::updateSurfaceStatusWithCurvature(
        stepper, state, *target, 0, Direction::Forward(),
        BoundaryTolerance::None(), s_onSurfaceTolerance,
        ConstrainedStep::Type::Navigator, getDummyLogger());
  };
  const auto field = std::make_shared<ConstantBField>(Vector3(0.005, 0., 0.));
  EigenStepper<> stepper(field);
  auto state = stepper.makeState(EigenStepper<>::Options(gctx, mfContext));
  stepper.initialize(state, start);
  const auto corrected = check(stepper, state);
  BOOST_REQUIRE(corrected.ok());
  BOOST_CHECK(*corrected == IntersectionStatus::unreachable);

  LegacySurfaceStatusStepper legacy(field);
  auto legacyState = legacy.makeState(AtlasStepper::Options(gctx, mfContext));
  legacy.initialize(legacyState, start);
  const auto original = check(legacy, legacyState);
  BOOST_REQUIRE(original.ok());
  BOOST_CHECK(*original == IntersectionStatus::unreachable);

  SurfaceStatusCurvatureRecorder recorder(field);
  auto recordingState =
      recorder.makeState(AtlasStepper::Options(gctx, mfContext));
  recorder.initialize(recordingState, start);
  const auto recorded = check(recorder, recordingState);
  BOOST_REQUIRE(recorded.ok());
  BOOST_CHECK(*recorded == IntersectionStatus::unreachable);
  BOOST_CHECK_SMALL(
      (recorder.receivedCurvature - Vector3(0., 0.005, 0.)).norm(), 1e-12);

  EigenStepper<> failing(std::make_shared<FailingCurvatureField>());
  auto failingState =
      failing.makeState(EigenStepper<>::Options(gctx, mfContext));
  failing.initialize(failingState, start);
  const auto error = check(failing, failingState);
  BOOST_REQUIRE(!error.ok());
  BOOST_CHECK(error.error() ==
              make_error_code(MagneticFieldError::OutOfBounds));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

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
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/detail/AlignmentHelper.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsAlignment/Geometry/AlignableStructure.hpp"
#include "ActsAlignment/Geometry/AlignmentHierarchy.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/detail/AlignmentEngine.hpp"
#include "ActsAlignment/Kernel/detail/AlignmentParameterDirectory.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <array>
#include <memory>
#include <stdexcept>
#include <variant>
#include <vector>

namespace {

using namespace Acts;
using namespace ActsAlignment;
using namespace Acts::UnitLiterals;

// Utility to create a dummy surface and detector element
auto makeDummyElement(const Acts::Transform3& trafo) {
  auto bounds = std::make_shared<const Acts::RectangleBounds>(10_mm, 10_mm);
  return std::make_shared<ActsTests::DetectorElementStub>(trafo, bounds, 1_um);
}

// Dummy options for the test
struct DummyFitOptions {
  Acts::GeometryContext geoContext =
      Acts::GeometryContext::dangerouslyDefaultConstruct();
  const Acts::Surface* referenceSurface = nullptr;
};

// Unique track type for hierarchy validation tests to avoid name collisions
struct HierarchyValidationTestTrack {
  std::size_t tipIndex() const { return 0; }
};

// Dummy fitter for alignment
struct DummyFitter {
  template <typename it_t, typename start_parameters_t, typename fit_options_t,
            typename track_container_t>
  Acts::Result<HierarchyValidationTestTrack> fit(
      it_t /*begin*/, it_t /*end*/,
      const start_parameters_t& /*startParameters*/,
      const fit_options_t& /*fitOptions*/,
      track_container_t& /*tracks*/) const {
    return Acts::Result<HierarchyValidationTestTrack>::success(
        HierarchyValidationTestTrack{});
  }
};

struct DerivativeFrames {
  Transform3 structure;
  Transform3 surface;
};

class ContextualDerivativeElement final
    : public ActsTests::DetectorElementStub {
 public:
  ContextualDerivativeElement()
      : DetectorElementStub(
            Transform3::Identity(),
            std::make_shared<const RectangleBounds>(1000_mm, 1000_mm), 1_um) {}

  const Transform3& localToGlobalTransform(
      const GeometryContext& gctx) const override {
    return gctx.get<DerivativeFrames>().surface;
  }
};

struct StructureDerivativeFixture {
  DerivativeFrames frames;
  ContextualDerivativeElement element;
  std::shared_ptr<AlignableStructure> structure;
  Vector3 trackOrigin = Vector3::Zero();
  Vector3 direction = Vector3::Zero();
  FreeVector pathDerivative = FreeVector::Zero();
  Matrix<2, eBoundSize> projector = Matrix<2, eBoundSize>::Zero();

  StructureDerivativeFixture() {
    frames.structure = Translation3(100_mm, -40_mm, 300_mm) *
                       AngleAxis3(0.6_rad, Vector3(1., 2., -1.).normalized());
    frames.surface = frames.structure * Translation3(30_mm, -20_mm, 50_mm) *
                     AngleAxis3(-0.8_rad, Vector3(2., -1., 3.).normalized());
    structure = std::make_shared<AlignableStructure>(
        GeometryIdentifier().withVolume(1), [](const GeometryContext& gctx) {
          return gctx.get<DerivativeFrames>().structure;
        });
    structure->addSurface(element.surface().getSharedPtr());
    structure->alignmentMask() = AlignmentMask::All;
    trackOrigin = frames.surface * Vector3(3_mm, -7_mm, 0.);
    direction =
        (frames.surface.rotation() * Vector3(0.3, -0.2, 1.)).normalized();
    pathDerivative.segment<3>(eFreePos0) = direction;
    // Reverse measurement order to exercise the bound-to-measurement
    // projection.
    projector(0, eBoundLoc1) = 1.;
    projector(1, eBoundLoc0) = 1.;
  }

  Vector3 intersection(const DerivativeFrames& state) const {
    const Vector3 normal = state.surface.rotation().col(2);
    const double path = normal.dot(state.surface.translation() - trackOrigin) /
                        normal.dot(direction);
    return trackOrigin + path * direction;
  }

  Vector2 residual(const DerivativeFrames& state) const {
    const GeometryContext gctx{state};
    BoundVector prediction = BoundVector::Zero();
    prediction.segment<2>(eBoundLoc0) =
        element.surface()
            .globalToLocal(gctx, intersection(state), direction)
            .value();
    return Vector2(1_mm, -2_mm) - projector * prediction;
  }

  DerivativeFrames moved(const DerivativeFrames& state, std::size_t parameter,
                         double step) const {
    Transform3 delta = Transform3::Identity();
    if (parameter < eAlignmentRotation0) {
      delta.translation()[parameter] = step;
    } else {
      delta.linear() =
          AngleAxis3(step, Vector3::Unit(parameter - eAlignmentRotation0))
              .toRotationMatrix();
    }
    const Transform3 current = state.structure * delta;
    return {current, current * state.structure.inverse() * state.surface};
  }
};

}  // namespace

BOOST_AUTO_TEST_CASE(StructureResidualDerivative) {
  const StructureDerivativeFixture fixture;
  const AlignmentHierarchy hierarchy({fixture.structure});
  BOOST_REQUIRE(hierarchy.validate().ok());
  const auto* owner = hierarchy.structureFor(fixture.element.surface());
  BOOST_REQUIRE_EQUAL(owner, fixture.structure.get());

  // Both contexts are displaced from the nominal element and have different
  // structure/member poses. Hold the same global straight track in both.
  const Transform3 motion =
      Translation3(-12_mm, 8_mm, 5_mm) * AngleAxis3(-0.2_rad, Vector3::UnitY());
  const std::array states{fixture.frames,
                          DerivativeFrames{motion * fixture.frames.structure,
                                           motion * fixture.frames.surface}};
  for (std::size_t contextIndex = 0; contextIndex < states.size();
       ++contextIndex) {
    const auto& state = states[contextIndex];
    const GeometryContext gctx{state};
    const Transform3 frame = owner->localToGlobalTransform(gctx);
    const Matrix<2, eAlignmentSize> analytical =
        -fixture.projector * ActsAlignment::detail::alignmentToBoundDerivative(
                                 gctx, fixture.element.surface(),
                                 fixture.intersection(state), fixture.direction,
                                 fixture.pathDerivative, owner->alignmentMask(),
                                 &frame);
    for (double scale : std::array{1e-5, 2e-6}) {
      for (std::size_t parameter = 0; parameter < eAlignmentSize; ++parameter) {
        BOOST_TEST_CONTEXT("context=" << contextIndex << ", parameter="
                                      << parameter << ", scale=" << scale) {
          const double step =
              scale * (parameter < eAlignmentRotation0 ? 1_mm : 1_rad);
          const Vector2 numerical =
              (fixture.residual(fixture.moved(state, parameter, step)) -
               fixture.residual(fixture.moved(state, parameter, -step))) /
              (2. * step);
          CHECK_CLOSE_ABS(analytical.col(parameter), numerical, 1e-6);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(StructureDerivativeMaskAfterConversion) {
  const StructureDerivativeFixture fixture;
  const GeometryContext gctx{fixture.frames};
  const Transform3 frame = fixture.structure->localToGlobalTransform(gctx);
  const auto& surface = fixture.element.surface();
  const Vector3 position = fixture.intersection(fixture.frames);
  const AlignmentMask mask = AlignmentMask::Center0 | AlignmentMask::Rotation2;
  const AlignmentToBoundMatrix masked =
      ActsAlignment::detail::alignmentToBoundDerivative(
          gctx, surface, position, fixture.direction, fixture.pathDerivative,
          mask, &frame);
  for (std::size_t parameter = 0; parameter < eAlignmentSize; ++parameter) {
    if (parameter == eAlignmentCenter0 || parameter == eAlignmentRotation2) {
      const double step =
          2e-6 * (parameter < eAlignmentRotation0 ? 1_mm : 1_rad);
      const Vector2 numerical =
          (fixture.residual(fixture.moved(fixture.frames, parameter, step)) -
           fixture.residual(fixture.moved(fixture.frames, parameter, -step))) /
          (2. * step);
      const Vector2 analytical = -fixture.projector * masked.col(parameter);
      CHECK_CLOSE_ABS(analytical, numerical, 1e-6);
    } else {
      BOOST_CHECK_EQUAL(masked.col(parameter).squaredNorm(), 0.);
    }
  }
  const AlignmentToBoundMatrix fixed =
      ActsAlignment::detail::alignmentToBoundDerivative(
          gctx, surface, position, fixture.direction, fixture.pathDerivative,
          AlignmentMask::None, &frame);
  BOOST_CHECK_EQUAL(fixed.squaredNorm(), 0.);

  // Verify that this rotated geometry distinguishes premature masking and
  // the obsolete conversion of surface translations to global coordinates.
  const AlignmentToBoundMatrix raw = surface.alignmentToBoundDerivative(
      gctx, position, fixture.direction, fixture.pathDerivative);
  AlignmentToBoundMatrix prematureMask = raw;
  ActsAlignment::detail::resetAlignmentDerivative(prematureMask, mask);
  const AlignmentMatrix composite =
      Acts::detail::compositeToComponentJacobian(frame, fixture.frames.surface);
  const Vector2 wronglyMasked =
      -fixture.projector * prematureMask * composite.col(eAlignmentCenter0);
  AlignmentMatrix obsoleteAdapter = AlignmentMatrix::Identity();
  obsoleteAdapter.block<3, 3>(eAlignmentCenter0, eAlignmentCenter0) =
      fixture.frames.surface.rotation();
  const Vector2 wronglyAdapted = -fixture.projector * raw * obsoleteAdapter *
                                 composite.col(eAlignmentCenter0);
  const Vector2 correct = -fixture.projector * masked.col(eAlignmentCenter0);
  BOOST_CHECK_GT((wronglyMasked - correct).norm(), 1e-2);
  BOOST_CHECK_GT((wronglyAdapted - correct).norm(), 1e-2);
}

BOOST_AUTO_TEST_CASE(StandaloneDerivativeCompatibility) {
  const StructureDerivativeFixture fixture;
  const GeometryContext gctx{fixture.frames};
  const auto& surface = fixture.element.surface();
  const Vector3 position = fixture.intersection(fixture.frames);
  for (AlignmentMask mask :
       {AlignmentMask::All, AlignmentMask::None,
        AlignmentMask::Center1 | AlignmentMask::Rotation0}) {
    AlignmentToBoundMatrix expected = surface.alignmentToBoundDerivative(
        gctx, position, fixture.direction, fixture.pathDerivative);
    ActsAlignment::detail::resetAlignmentDerivative(expected, mask);
    const AlignmentToBoundMatrix actual =
        ActsAlignment::detail::alignmentToBoundDerivative(
            gctx, surface, position, fixture.direction, fixture.pathDerivative,
            mask);
    BOOST_CHECK_EQUAL((actual - expected).squaredNorm(), 0.);
  }
}

BOOST_AUTO_TEST_CASE(StructureWithoutFrame) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  const AlignableStructure structure(GeometryIdentifier().withVolume(1));
  BOOST_CHECK(!structure.hasTransform());
  BOOST_CHECK_THROW(structure.localToGlobalTransform(gctx), std::logic_error);
  BOOST_CHECK_THROW(AlignableStructure(GeometryIdentifier().withVolume(1),
                                       AlignableStructure::TransformProvider{}),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(FixedStructureFrame) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  const Transform3 expected = Translation3(100_mm, -30_mm, 200_mm) *
                              AngleAxis3(0.7_rad, Vector3::UnitY());
  Transform3 supplied = expected;
  auto parent = std::make_shared<AlignableStructure>(
      GeometryIdentifier().withVolume(1), supplied);
  supplied.setIdentity();
  BOOST_CHECK(parent->hasTransform());
  CHECK_CLOSE_ABS(parent->localToGlobalTransform(gctx).matrix(),
                  expected.matrix(), 1e-14);

  // An explicitly global child frame is not composed with its parent's frame.
  const Transform3 childFrame = Translation3(110_mm, -20_mm, 220_mm) *
                                AngleAxis3(-0.3_rad, Vector3::UnitX());
  auto child = std::make_shared<AlignableStructure>(
      GeometryIdentifier().withVolume(2), childFrame);
  parent->addChild(child);
  auto element = makeDummyElement(Transform3::Identity());
  child->addSurface(element->surface().getSharedPtr());
  const AlignmentHierarchy hierarchy({parent});
  BOOST_REQUIRE(hierarchy.validate().ok());
  BOOST_REQUIRE_EQUAL(hierarchy.structureFor(element->surface()), child.get());
  CHECK_CLOSE_ABS(child->localToGlobalTransform(gctx).matrix(),
                  childFrame.matrix(), 1e-14);
  // Declaring or querying a frame does not move the nominal member surface.
  CHECK_CLOSE_ABS(element->surface().localToGlobalTransform(gctx).matrix(),
                  Transform3::Identity().matrix(), 1e-14);
}

BOOST_AUTO_TEST_CASE(ContextualStructureFrame) {
  struct FrameContext {
    Transform3 transform;
  };
  const AlignableStructure structure(
      GeometryIdentifier().withVolume(1), [](const GeometryContext& gctx) {
        return gctx.get<FrameContext>().transform;
      });
  const Transform3 first =
      Translation3(12_mm, -5_mm, 70_mm) * AngleAxis3(0.4_rad, Vector3::UnitZ());
  const Transform3 second = Translation3(-18_mm, 9_mm, 80_mm) *
                            AngleAxis3(-0.2_rad, Vector3::UnitY());
  const GeometryContext firstContext{FrameContext{first}};
  const GeometryContext secondContext{FrameContext{second}};
  BOOST_CHECK(structure.hasTransform());
  auto returned = structure.localToGlobalTransform(firstContext);
  CHECK_CLOSE_ABS(returned.matrix(), first.matrix(), 1e-14);
  returned.setIdentity();
  CHECK_CLOSE_ABS(structure.localToGlobalTransform(secondContext).matrix(),
                  second.matrix(), 1e-14);
  // A query in another context, or mutation of a returned value, must not
  // replace or modify the frame seen by the first context.
  CHECK_CLOSE_ABS(structure.localToGlobalTransform(firstContext).matrix(),
                  first.matrix(), 1e-14);
}

BOOST_AUTO_TEST_CASE(StructureFrameReplacement) {
  struct FrameContext {
    Transform3 transform;
  };
  const Transform3 first =
      Translation3(12_mm, -5_mm, 70_mm) * AngleAxis3(0.4_rad, Vector3::UnitZ());
  const Transform3 second = Translation3(-18_mm, 9_mm, 80_mm) *
                            AngleAxis3(-0.2_rad, Vector3::UnitY());
  const GeometryContext firstContext{FrameContext{first}};
  const GeometryContext secondContext{FrameContext{second}};
  auto structure =
      std::make_shared<AlignableStructure>(GeometryIdentifier().withVolume(1));
  auto element = makeDummyElement(Transform3::Identity());
  structure->addSurface(element->surface().getSharedPtr());
  const AlignmentHierarchy hierarchy({structure});
  BOOST_REQUIRE(hierarchy.validate().ok());

  BOOST_CHECK_THROW(
      structure->setTransform(AlignableStructure::TransformProvider{}),
      std::invalid_argument);
  BOOST_CHECK(!structure->hasTransform());

  Transform3 supplied = first;
  structure->setTransform(supplied);
  supplied.setIdentity();
  BOOST_REQUIRE(structure->hasTransform());
  CHECK_CLOSE_ABS(structure->localToGlobalTransform(firstContext).matrix(),
                  first.matrix(), 1e-14);

  const auto provider = [](const GeometryContext& gctx) {
    return gctx.get<FrameContext>().transform;
  };
  structure->setTransform(provider);
  CHECK_CLOSE_ABS(structure->localToGlobalTransform(firstContext).matrix(),
                  first.matrix(), 1e-14);
  CHECK_CLOSE_ABS(structure->localToGlobalTransform(secondContext).matrix(),
                  second.matrix(), 1e-14);
  BOOST_CHECK_THROW(
      structure->setTransform(AlignableStructure::TransformProvider{}),
      std::invalid_argument);
  CHECK_CLOSE_ABS(structure->localToGlobalTransform(secondContext).matrix(),
                  second.matrix(), 1e-14);

  structure->setTransform(first);
  CHECK_CLOSE_ABS(structure->localToGlobalTransform(secondContext).matrix(),
                  first.matrix(), 1e-14);
  BOOST_CHECK_THROW(
      structure->setTransform(AlignableStructure::TransformProvider{}),
      std::invalid_argument);
  CHECK_CLOSE_ABS(structure->localToGlobalTransform(secondContext).matrix(),
                  first.matrix(), 1e-14);

  // Frames can be attached after building the hierarchy without changing its
  // membership lookup or moving a member surface.
  BOOST_CHECK_EQUAL(hierarchy.structureFor(element->surface()),
                    structure.get());
  CHECK_CLOSE_ABS(
      element->surface().localToGlobalTransform(firstContext).matrix(),
      Transform3::Identity().matrix(), 1e-14);

  AlignableStructure lateContextual(GeometryIdentifier().withVolume(2));
  lateContextual.setTransform(provider);
  BOOST_CHECK(lateContextual.hasTransform());
  CHECK_CLOSE_ABS(lateContextual.localToGlobalTransform(secondContext).matrix(),
                  second.matrix(), 1e-14);
}

BOOST_AUTO_TEST_CASE(AlignmentRejectsOverlappingStructures) {
  auto element = makeDummyElement(Acts::Transform3::Identity());
  auto struct1 = std::make_shared<ActsAlignment::AlignableStructure>(
      Acts::GeometryIdentifier().withVolume(1));
  auto struct2 = std::make_shared<ActsAlignment::AlignableStructure>(
      Acts::GeometryIdentifier().withVolume(2));
  struct1->addSurface(element->surface().getSharedPtr());
  struct2->addSurface(element->surface().getSharedPtr());

  // Suppress the expected error log for this negative test only. This works
  // with both disabled and enabled log-failure thresholds, without changing
  // the global threshold used by other tests.
  ActsAlignment::Alignment alignEngine(DummyFitter{},
                                       Acts::getDummyLogger().clone());

  // Validation must reject the hierarchy even with no solver iterations.
  // Empty trajectories are never passed to the numerical solve.
  ActsAlignment::AlignmentOptions<DummyFitOptions> options(
      DummyFitOptions{}, {}, {element.get()}, 0.5, {5, 0.01}, 0, {},
      {struct1, struct2});
  const std::vector<std::vector<Acts::SourceLink>> trajs;
  const std::vector<Acts::BoundTrackParameters> params;

  const auto result = alignEngine.align(trajs, params, options);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK_EQUAL(result.error(),
                    ActsAlignment::AlignmentError::HierarchyValidationFailure);
}

BOOST_AUTO_TEST_CASE(AlignmentHierarchyHelper) {
  // Direct test of AlignmentHierarchy without standing up Alignment<Fitter>.
  auto el1 = makeDummyElement(Acts::Translation3(0, 0, 10_mm) *
                              Acts::Transform3::Identity());
  auto el2 = makeDummyElement(Acts::Translation3(0, 0, 20_mm) *
                              Acts::Transform3::Identity());
  auto el3 = makeDummyElement(Acts::Translation3(0, 0, 30_mm) *
                              Acts::Transform3::Identity());

  auto structA = std::make_shared<ActsAlignment::AlignableStructure>(
      Acts::GeometryIdentifier().withVolume(10));
  auto structB = std::make_shared<ActsAlignment::AlignableStructure>(
      Acts::GeometryIdentifier().withVolume(20));

  // structA owns el1 + el2; el3 floats as a standalone module.
  structA->addSurface(el1->surface().getSharedPtr());
  structA->addSurface(el2->surface().getSharedPtr());

  // --- Clean hierarchy ---
  {
    ActsAlignment::AlignmentHierarchy hierarchy({structA, structB});
    BOOST_CHECK(hierarchy.validate().ok());
    BOOST_CHECK_EQUAL(hierarchy.structureFor(*el1), structA.get());
    BOOST_CHECK_EQUAL(hierarchy.structureFor(*el2), structA.get());
    // el3 is a standalone floating module
    BOOST_CHECK_EQUAL(hierarchy.structureFor(*el3), nullptr);
  }

  // --- Surface-level lookup ---
  {
    ActsAlignment::AlignmentHierarchy hierarchy({structA, structB});
    BOOST_CHECK_EQUAL(hierarchy.structureFor(el1->surface()), structA.get());
    BOOST_CHECK_EQUAL(hierarchy.structureFor(el2->surface()), structA.get());
    BOOST_CHECK_EQUAL(hierarchy.structureFor(el3->surface()), nullptr);
  }

  // --- Overlap detected ---
  structB->addSurface(el1->surface().getSharedPtr());
  {
    ActsAlignment::AlignmentHierarchy hierarchy({structA, structB});
    const auto result = hierarchy.validate();
    BOOST_CHECK(!result.ok());
    BOOST_REQUIRE_EQUAL(result.overlapping.size(), 1u);
    BOOST_CHECK_EQUAL(result.overlapping.front(), &el1->surface());
  }

  // --- Empty hierarchy ---
  {
    ActsAlignment::AlignmentHierarchy hierarchy({});
    BOOST_CHECK(hierarchy.validate().ok());
    BOOST_CHECK_EQUAL(hierarchy.structureFor(*el1), nullptr);
  }
}

BOOST_AUTO_TEST_CASE(NestedHierarchyRegistration) {
  auto el1 = makeDummyElement(Acts::Translation3(0, 0, 10_mm) *
                              Acts::Transform3::Identity());
  auto el2 = makeDummyElement(Acts::Translation3(0, 0, 20_mm) *
                              Acts::Transform3::Identity());
  auto el3 = makeDummyElement(Acts::Translation3(0, 0, 30_mm) *
                              Acts::Transform3::Identity());
  auto parent = std::make_shared<ActsAlignment::AlignableStructure>(
      Acts::GeometryIdentifier().withVolume(1));
  auto child = std::make_shared<ActsAlignment::AlignableStructure>(
      Acts::GeometryIdentifier().withVolume(2));
  auto grandchild = std::make_shared<ActsAlignment::AlignableStructure>(
      Acts::GeometryIdentifier().withVolume(3));
  parent->addSurface(el1->surface().getSharedPtr());
  child->addSurface(el2->surface().getSharedPtr());
  grandchild->addSurface(el3->surface().getSharedPtr());
  parent->addChild(child);
  child->addChild(grandchild);

  // Register only the root. Descendants must be visited recursively, and
  // each surface must retain its immediate owner instead of the root.
  {
    ActsAlignment::AlignmentHierarchy hierarchy({parent});
    BOOST_CHECK(hierarchy.validate().ok());
    BOOST_CHECK_EQUAL(hierarchy.structureFor(el1->surface()), parent.get());
    BOOST_CHECK_EQUAL(hierarchy.structureFor(el2->surface()), child.get());
    BOOST_CHECK_EQUAL(hierarchy.structureFor(el3->surface()), grandchild.get());
  }

  // Duplicate ownership must also be detected across hierarchy levels.
  grandchild->addSurface(el1->surface().getSharedPtr());
  {
    ActsAlignment::AlignmentHierarchy hierarchy({parent});
    const auto result = hierarchy.validate();
    BOOST_CHECK(!result.ok());
    BOOST_REQUIRE_EQUAL(result.overlapping.size(), 1u);
    BOOST_CHECK_EQUAL(result.overlapping.front(), &el1->surface());
  }
}

// TODO: re-enable once detectMaskConflicts() is implemented. The function is
// currently stubbed to return empty (see AlignmentHierarchy.hpp) because a
// correct conflict check requires the chain-rule Jacobian to compare DoFs
// across coordinate frames.
BOOST_AUTO_TEST_CASE(MaskConflictDetection, *boost::unit_test::disabled()) {
  using ActsAlignment::AlignmentMask;

  auto el1 = makeDummyElement(Acts::Translation3(0, 0, 10_mm) *
                              Acts::Transform3::Identity());
  auto el2 = makeDummyElement(Acts::Translation3(0, 0, 20_mm) *
                              Acts::Transform3::Identity());
  auto el3 = makeDummyElement(Acts::Translation3(0, 0, 30_mm) *
                              Acts::Transform3::Identity());

  auto parent = std::make_shared<ActsAlignment::AlignableStructure>(
      Acts::GeometryIdentifier().withVolume(1));
  parent->addSurface(el1->surface().getSharedPtr());
  parent->addSurface(el2->surface().getSharedPtr());

  std::vector<Acts::SurfacePlacementBase*> floatingModules = {
      el1.get(), el2.get(), el3.get()};

  ActsAlignment::AlignmentHierarchy hierarchy({parent});

  // Parent floats Center0. Module iteration mask = All → conflict on el1, el2.
  parent->alignmentMask() = AlignmentMask::Center0;
  {
    const auto conflicts =
        hierarchy.detectMaskConflicts(AlignmentMask::All, floatingModules);
    BOOST_CHECK_EQUAL(conflicts.size(), 2u);
    for (const auto& c : conflicts) {
      BOOST_CHECK_EQUAL(c.structure, parent.get());
      BOOST_CHECK(c.conflictingBits == AlignmentMask::Center0);
    }
  }

  // Disjoint masks: parent = Center0, module mask = Rotation0 → no conflict.
  {
    const auto conflicts = hierarchy.detectMaskConflicts(
        AlignmentMask::Rotation0, floatingModules);
    BOOST_CHECK(conflicts.empty());
  }

  // Parent with None mask → no conflict regardless of module mask.
  parent->alignmentMask() = AlignmentMask::None;
  {
    const auto conflicts =
        hierarchy.detectMaskConflicts(AlignmentMask::All, floatingModules);
    BOOST_CHECK(conflicts.empty());
  }

  // Standalone module (el3) is never in a structure → never a conflict.
  parent->alignmentMask() = AlignmentMask::All;
  {
    const auto conflicts =
        hierarchy.detectMaskConflicts(AlignmentMask::All, {el3.get()});
    BOOST_CHECK(conflicts.empty());
  }
}

namespace {
namespace AlignDetail = ActsAlignment::detail;

auto makeParameterDirectory(
    const std::vector<std::shared_ptr<AlignableStructure>>& registered,
    const std::vector<std::shared_ptr<AlignableStructure>>& selected,
    const std::vector<SurfacePlacementBase*>& modules = {},
    AlignmentMask moduleMask = AlignmentMask::All) {
  return AlignDetail::makeAlignmentParameterDirectory(registered, selected,
                                                      modules, moduleMask);
}

void checkSelectionError(
    const Acts::Result<AlignDetail::AlignmentParameterDirectory,
                       AlignDetail::AlignmentSelectionError>& result,
    AlignmentError code, const std::string& message) {
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error().code == code);
  BOOST_CHECK_NE(result.error().message.find(message), std::string::npos);
}

struct ParameterDirectoryFixture {
  std::shared_ptr<ActsTests::DetectorElementStub> first =
      makeDummyElement(Transform3::Identity());
  std::shared_ptr<ActsTests::DetectorElementStub> second =
      makeDummyElement(Transform3::Identity());
  std::shared_ptr<ActsTests::DetectorElementStub> module =
      makeDummyElement(Transform3::Identity());
  std::shared_ptr<ActsTests::DetectorElementStub> other =
      makeDummyElement(Transform3::Identity());
  // Reverse the IDs to check that numbering follows selection order.
  std::shared_ptr<AlignableStructure> group =
      std::make_shared<AlignableStructure>(GeometryIdentifier().withVolume(2));
  std::shared_ptr<AlignableStructure> otherGroup =
      std::make_shared<AlignableStructure>(GeometryIdentifier().withVolume(1));

  ParameterDirectoryFixture() {
    group->alignmentMask() = AlignmentMask::All;
    group->addSurface(first->surface().getSharedPtr());
    group->addSurface(second->surface().getSharedPtr());
    otherGroup->alignmentMask() = AlignmentMask::Center2;
    otherGroup->addSurface(other->surface().getSharedPtr());
  }
};
}  // namespace

BOOST_FIXTURE_TEST_CASE(DirectoryRegistrationPreservesAssignments,
                        ParameterDirectoryFixture) {
  AlignDetail::AlignmentParameterDirectory directory;
  const auto modulePosition =
      directory.addAlignable({module.get(), {}}, AlignmentMask::None);
  const auto structurePosition = directory.addAlignable(
      {group.get(), {}}, AlignmentMask::Center1 | AlignmentMask::Rotation2);
  BOOST_CHECK_EQUAL(modulePosition, 0u);
  BOOST_CHECK_EQUAL(structurePosition, 1u);
  BOOST_REQUIRE_EQUAL(directory.alignables.size(), 2u);
  BOOST_CHECK_EQUAL(directory.numParameters, 2u);
  for (const auto& parameter :
       directory.alignables[modulePosition].parameterIndices) {
    BOOST_CHECK(!parameter);
  }
  const auto& indices =
      directory.alignables[structurePosition].parameterIndices;
  BOOST_REQUIRE(indices[eAlignmentCenter1]);
  BOOST_REQUIRE(indices[eAlignmentRotation2]);
  BOOST_CHECK_EQUAL(*indices[eAlignmentCenter1], 0u);
  BOOST_CHECK_EQUAL(*indices[eAlignmentRotation2], 1u);
  BOOST_CHECK(!indices[eAlignmentCenter0]);

  BOOST_CHECK(directory.registerPlacement(module.get(), modulePosition));
  BOOST_CHECK(directory.registerPlacement(first.get(), structurePosition));
  BOOST_CHECK(directory.registerPlacement(second.get(), structurePosition));
  // Reject repeated registration without replacing the original association.
  BOOST_CHECK(!directory.registerPlacement(module.get(), structurePosition));
  BOOST_CHECK(!directory.registerPlacement(first.get(), structurePosition));
  BOOST_CHECK_EQUAL(directory.placementToAlignable.size(), 3u);
  BOOST_CHECK_EQUAL(directory.placementToAlignable.at(module.get()),
                    modulePosition);
  BOOST_CHECK_EQUAL(directory.placementToAlignable.at(first.get()),
                    structurePosition);
  BOOST_CHECK_EQUAL(directory.placementToAlignable.at(second.get()),
                    structurePosition);
  BOOST_CHECK_EQUAL(directory.numParameters, 2u);
}

BOOST_FIXTURE_TEST_CASE(DirectoryRenumbersCopiedUnits,
                        ParameterDirectoryFixture) {
  const auto original = makeParameterDirectory({group}, {group});
  BOOST_REQUIRE(original.ok());
  AlignDetail::AlignmentParameterDirectory directory;
  directory.addAlignable({module.get(), {}}, AlignmentMask::Center0);
  const auto position = directory.addAlignable(original->alignables.front(),
                                               AlignmentMask::Rotation1);
  BOOST_REQUIRE_EQUAL(directory.alignables.size(), 2u);
  BOOST_CHECK_EQUAL(position, 1u);
  BOOST_CHECK_EQUAL(directory.numParameters, 2u);
  const auto& indices = directory.alignables[position].parameterIndices;
  BOOST_REQUIRE(indices[eAlignmentRotation1]);
  BOOST_CHECK_EQUAL(*indices[eAlignmentRotation1], 1u);
  for (std::size_t dof = 0; dof < eAlignmentSize; ++dof) {
    if (dof != eAlignmentRotation1) {
      BOOST_CHECK(!indices[dof]);
    }
    BOOST_REQUIRE(original->alignables.front().parameterIndices[dof]);
    BOOST_CHECK_EQUAL(*original->alignables.front().parameterIndices[dof], dof);
  }
}

BOOST_FIXTURE_TEST_CASE(MembersShareParameterIndices,
                        ParameterDirectoryFixture) {
  const auto directory = makeParameterDirectory({group}, {group});
  BOOST_REQUIRE(directory.ok());
  BOOST_REQUIRE_EQUAL(directory->alignables.size(), 1u);
  BOOST_CHECK_EQUAL(directory->numParameters, 6u);
  BOOST_CHECK_EQUAL(directory->placementToAlignable.size(), 2u);
  BOOST_CHECK_EQUAL(directory->placementToAlignable.at(first.get()), 0u);
  BOOST_CHECK_EQUAL(directory->placementToAlignable.at(second.get()), 0u);
  const auto& entry = directory->alignables.front();
  BOOST_CHECK(std::get<const AlignableStructure*>(entry.alignable) ==
              group.get());
  for (std::size_t dof = 0; dof < eAlignmentSize; ++dof) {
    BOOST_REQUIRE(entry.parameterIndices[dof].has_value());
    BOOST_CHECK_EQUAL(*entry.parameterIndices[dof], dof);
  }
  BOOST_CHECK_EQUAL(directory->placementToAlignable.at(first.get()),
                    directory->placementToAlignable.at(second.get()));
}

BOOST_FIXTURE_TEST_CASE(MixedSelectionAndParameterOrder,
                        ParameterDirectoryFixture) {
  group->alignmentMask() = AlignmentMask::Center0 | AlignmentMask::Rotation2;
  const auto directory = makeParameterDirectory(
      {otherGroup, group}, {group, otherGroup}, {module.get()},
      AlignmentMask::Center1 | AlignmentMask::Rotation0);
  BOOST_REQUIRE(directory.ok());
  BOOST_REQUIRE_EQUAL(directory->alignables.size(), 3u);
  BOOST_CHECK_EQUAL(directory->numParameters, 5u);
  BOOST_CHECK(std::get<const SurfacePlacementBase*>(
                  directory->alignables[0].alignable) == module.get());
  BOOST_CHECK_EQUAL(directory->placementToAlignable.at(module.get()), 0u);
  BOOST_CHECK_EQUAL(directory->placementToAlignable.at(first.get()), 1u);
  BOOST_CHECK_EQUAL(directory->placementToAlignable.at(other.get()), 2u);
  const auto& moduleIndices = directory->alignables[0].parameterIndices;
  const auto& structureIndices = directory->alignables[1].parameterIndices;
  BOOST_REQUIRE(moduleIndices[eAlignmentCenter1]);
  BOOST_REQUIRE(moduleIndices[eAlignmentRotation0]);
  BOOST_REQUIRE(structureIndices[eAlignmentCenter0]);
  BOOST_REQUIRE(structureIndices[eAlignmentRotation2]);
  BOOST_REQUIRE(directory->alignables[2].parameterIndices[eAlignmentCenter2]);
  BOOST_CHECK_EQUAL(*moduleIndices[eAlignmentCenter1], 0u);
  BOOST_CHECK_EQUAL(*moduleIndices[eAlignmentRotation0], 1u);
  BOOST_CHECK_EQUAL(*structureIndices[eAlignmentCenter0], 2u);
  BOOST_CHECK_EQUAL(*structureIndices[eAlignmentRotation2], 3u);
  BOOST_CHECK_EQUAL(
      *directory->alignables[2].parameterIndices[eAlignmentCenter2], 4u);
  BOOST_CHECK(!moduleIndices[eAlignmentCenter0]);
  BOOST_CHECK(!structureIndices[eAlignmentCenter1]);

  const auto reversed =
      makeParameterDirectory({group, otherGroup}, {otherGroup, group});
  BOOST_REQUIRE(reversed.ok());
  BOOST_CHECK_EQUAL(reversed->placementToAlignable.at(other.get()), 0u);
  BOOST_CHECK_EQUAL(reversed->placementToAlignable.at(first.get()), 1u);
  BOOST_CHECK_EQUAL(
      *reversed->alignables[0].parameterIndices[eAlignmentCenter2], 0u);
  BOOST_CHECK_EQUAL(
      *reversed->alignables[1].parameterIndices[eAlignmentCenter0], 1u);
}

BOOST_FIXTURE_TEST_CASE(InactiveSelectionsAndEffectiveMasks,
                        ParameterDirectoryFixture) {
  auto metadata = makeParameterDirectory({group}, {});
  BOOST_REQUIRE(metadata.ok());
  BOOST_CHECK(metadata->alignables.empty());
  BOOST_CHECK_EQUAL(metadata->numParameters, 0u);

  group->alignmentMask() = AlignmentMask::None;
  auto inactive = makeParameterDirectory({group}, {group}, {first.get()},
                                         AlignmentMask::None);
  BOOST_REQUIRE(inactive.ok());
  BOOST_REQUIRE_EQUAL(inactive->alignables.size(), 1u);
  BOOST_CHECK_EQUAL(inactive->numParameters, 0u);
  BOOST_CHECK_EQUAL(inactive->placementToAlignable.at(first.get()), 0u);
  BOOST_CHECK(!inactive->placementToAlignable.contains(second.get()));
  for (const auto& parameter : inactive->alignables[0].parameterIndices) {
    BOOST_CHECK(!parameter);
  }
  // Changing the mask requires a new directory.
  auto enabled = makeParameterDirectory({group}, {group}, {first.get()},
                                        AlignmentMask::Rotation1);
  BOOST_REQUIRE(enabled.ok());
  BOOST_REQUIRE(enabled->alignables[0].parameterIndices[eAlignmentRotation1]);
  BOOST_CHECK_EQUAL(
      *enabled->alignables[0].parameterIndices[eAlignmentRotation1], 0u);
  BOOST_CHECK_EQUAL(enabled->numParameters, 1u);
  BOOST_CHECK_EQUAL(inactive->numParameters, 0u);

  auto empty = makeParameterDirectory({}, {});
  BOOST_REQUIRE(empty.ok());
  BOOST_CHECK(empty->placementToAlignable.empty());
  BOOST_CHECK_EQUAL(empty->numParameters, 0u);
}

BOOST_FIXTURE_TEST_CASE(StandaloneParameterOrder, ParameterDirectoryFixture) {
  // An unselected structure leaves its members available as standalone modules.
  const auto directory = makeParameterDirectory(
      {group}, {}, {second.get(), first.get()}, AlignmentMask::Rotation2);
  BOOST_REQUIRE(directory.ok());
  BOOST_REQUIRE_EQUAL(directory->alignables.size(), 2u);
  BOOST_CHECK_EQUAL(directory->numParameters, 2u);
  BOOST_CHECK_EQUAL(directory->placementToAlignable.at(second.get()), 0u);
  BOOST_CHECK_EQUAL(directory->placementToAlignable.at(first.get()), 1u);
  BOOST_REQUIRE(directory->alignables[0].parameterIndices[eAlignmentRotation2]);
  BOOST_REQUIRE(directory->alignables[1].parameterIndices[eAlignmentRotation2]);
  BOOST_CHECK_EQUAL(
      *directory->alignables[0].parameterIndices[eAlignmentRotation2], 0u);
  BOOST_CHECK_EQUAL(
      *directory->alignables[1].parameterIndices[eAlignmentRotation2], 1u);
}

BOOST_FIXTURE_TEST_CASE(SelectionDoesNotResolveGeometry,
                        ParameterDirectoryFixture) {
  BOOST_CHECK(!group->hasTransform());
  auto withoutFrame = makeParameterDirectory({group}, {group});
  BOOST_REQUIRE(withoutFrame.ok());
  group->setTransform([](const GeometryContext&) -> Transform3 {
    throw std::runtime_error("parameter numbering must not resolve a frame");
  });
  // Numbering depends on the mask, not on numerical constraints.
  group->constraints()[eAlignmentCenter0] = 0.;
  auto directory = makeParameterDirectory({group}, {group});
  BOOST_REQUIRE(directory.ok());
  BOOST_CHECK_EQUAL(directory->numParameters, 6u);
}

BOOST_FIXTURE_TEST_CASE(ParameterDirectoryDoesNotOwnAlignables,
                        ParameterDirectoryFixture) {
  const auto structureOwners = group.use_count();
  const auto moduleOwners = module.use_count();
  const auto directory =
      makeParameterDirectory({group}, {group}, {module.get()});
  BOOST_REQUIRE(directory.ok());
  BOOST_CHECK_EQUAL(group.use_count(), structureOwners);
  BOOST_CHECK_EQUAL(module.use_count(), moduleOwners);
}

BOOST_FIXTURE_TEST_CASE(SelectedLeavesAndInactiveParents,
                        ParameterDirectoryFixture) {
  auto parent = std::make_shared<AlignableStructure>(GeometryIdentifier());
  parent->addChild(group);
  auto directory = makeParameterDirectory({parent}, {parent, group});
  BOOST_REQUIRE(directory.ok());
  BOOST_CHECK_EQUAL(directory->alignables.size(), 1u);
  parent->alignmentMask() = AlignmentMask::All;
  checkSelectionError(makeParameterDirectory({parent}, {parent}),
                      AlignmentError::InvalidAlignmentSelection, "leaf");
  checkSelectionError(makeParameterDirectory({parent}, {parent, group}),
                      AlignmentError::InvalidAlignmentSelection, "leaf");
}

BOOST_FIXTURE_TEST_CASE(InvalidParameterSelections, ParameterDirectoryFixture) {
  checkSelectionError(makeParameterDirectory({group}, {nullptr}),
                      AlignmentError::InvalidAlignmentSelection, "null");
  checkSelectionError(makeParameterDirectory({group}, {otherGroup}),
                      AlignmentError::InvalidAlignmentSelection,
                      "unregistered");
  checkSelectionError(makeParameterDirectory({group}, {group, group}),
                      AlignmentError::InvalidAlignmentSelection, "Duplicate");
  checkSelectionError(makeParameterDirectory({}, {}, {nullptr}),
                      AlignmentError::InvalidAlignmentSelection, "Null");
  checkSelectionError(
      makeParameterDirectory({}, {}, {module.get(), module.get()}),
      AlignmentError::InvalidAlignmentSelection, "Duplicate");
  // A structure member cannot also be selected separately, regardless of masks.
  group->alignmentMask() = AlignmentMask::Center0;
  for (const auto mask : {AlignmentMask::Rotation0, AlignmentMask::None}) {
    checkSelectionError(
        makeParameterDirectory({group}, {group}, {first.get()}, mask),
        AlignmentError::InvalidAlignmentSelection, "Overlapping");
  }
  group->alignmentMask() = AlignmentMask::None;
  checkSelectionError(makeParameterDirectory({group}, {group, group}),
                      AlignmentError::InvalidAlignmentSelection, "Duplicate");
  auto empty = std::make_shared<AlignableStructure>(GeometryIdentifier());
  empty->alignmentMask() = AlignmentMask::All;
  checkSelectionError(makeParameterDirectory({empty}, {empty}),
                      AlignmentError::InvalidAlignmentSelection, "nonempty");
  empty->addSurface(Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(1_mm, 1_mm)));
  // A placement is required only when the structure is selected and active.
  const auto unselected = makeParameterDirectory({empty}, {});
  BOOST_REQUIRE(unselected.ok());
  BOOST_CHECK(unselected->alignables.empty());
  BOOST_CHECK_EQUAL(unselected->numParameters, 0u);
  checkSelectionError(makeParameterDirectory({empty}, {empty}),
                      AlignmentError::InvalidAlignmentSelection, "placement");
  empty->alignmentMask() = AlignmentMask::None;
  const auto inactive = makeParameterDirectory({empty}, {empty});
  BOOST_REQUIRE(inactive.ok());
  BOOST_CHECK(inactive->alignables.empty());
  BOOST_CHECK_EQUAL(inactive->numParameters, 0u);
}

BOOST_FIXTURE_TEST_CASE(InvalidHierarchyForParameterDirectory,
                        ParameterDirectoryFixture) {
  checkSelectionError(makeParameterDirectory({nullptr}, {}),
                      AlignmentError::HierarchyValidationFailure, "Null");
  checkSelectionError(makeParameterDirectory({group, group}, {}),
                      AlignmentError::HierarchyValidationFailure, "repeated");
  auto parent = std::make_shared<AlignableStructure>(GeometryIdentifier());
  parent->addChild(group);
  parent->addChild(group);
  checkSelectionError(makeParameterDirectory({parent}, {}),
                      AlignmentError::HierarchyValidationFailure, "repeated");
  auto invalid = std::make_shared<AlignableStructure>(GeometryIdentifier());
  invalid->addChild(nullptr);
  checkSelectionError(makeParameterDirectory({invalid}, {}),
                      AlignmentError::HierarchyValidationFailure, "Null");
  invalid = std::make_shared<AlignableStructure>(GeometryIdentifier());
  invalid->addSurface(nullptr);
  checkSelectionError(makeParameterDirectory({invalid}, {}),
                      AlignmentError::HierarchyValidationFailure, "Null");
  otherGroup->addSurface(first->surface().getSharedPtr());
  checkSelectionError(makeParameterDirectory({group, otherGroup}, {}),
                      AlignmentError::HierarchyValidationFailure, "placement");
  group->addSurface(first->surface().getSharedPtr());
  checkSelectionError(makeParameterDirectory({group}, {}),
                      AlignmentError::HierarchyValidationFailure, "placement");

  // Non-owning aliases let us test a cycle without creating an ownership leak.
  AlignableStructure cycleA{GeometryIdentifier()}, cycleB{GeometryIdentifier()};
  auto a = std::shared_ptr<AlignableStructure>(&cycleA, [](auto*) {});
  auto b = std::shared_ptr<AlignableStructure>(&cycleB, [](auto*) {});
  cycleA.addChild(b);
  cycleB.addChild(a);
  checkSelectionError(makeParameterDirectory({a}, {}),
                      AlignmentError::HierarchyValidationFailure, "Cycle");
}

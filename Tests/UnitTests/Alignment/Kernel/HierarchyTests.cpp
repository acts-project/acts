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
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsAlignment/Geometry/AlignableStructure.hpp"
#include "ActsAlignment/Geometry/AlignmentHierarchy.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/detail/AlignmentTarget.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
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

}  // namespace

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

const AlignDetail::AlignmentUpdateSupport updateSupported =
    [](const SurfacePlacementBase&) { return true; };

auto validateSelection(
    const std::vector<std::shared_ptr<AlignableStructure>>& registered,
    const std::vector<std::shared_ptr<AlignableStructure>>& active) {
  return AlignDetail::validateAlignmentSelection(registered, active);
}

auto prepareTargets(
    const GeometryContext& gctx,
    const std::vector<std::shared_ptr<AlignableStructure>>& registered,
    const std::vector<std::shared_ptr<AlignableStructure>>& active,
    const std::vector<SurfacePlacementBase*>& modules = {},
    const AlignDetail::AlignmentUpdateSupport& support = updateSupported) {
  return AlignDetail::prepareAlignmentTargets(gctx, registered, active, modules,
                                              support);
}

template <typename T>
void checkTargetError(
    const Acts::Result<T, AlignDetail::AlignmentTargetError>& r,
    AlignmentError code, const std::string& diagnostic) {
  BOOST_REQUIRE(!r.ok());
  BOOST_CHECK(r.error().code == code);
  BOOST_CHECK_NE(r.error().message.find(diagnostic), std::string::npos);
}

struct TargetFixture {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  std::shared_ptr<ActsTests::DetectorElementStub> first =
      makeDummyElement(Transform3::Identity());
  std::shared_ptr<ActsTests::DetectorElementStub> second = makeDummyElement(
      Translation3(10_mm, 20_mm, 30_mm) * AngleAxis3(0.3, Vector3::UnitY()));
  std::shared_ptr<ActsTests::DetectorElementStub> module =
      makeDummyElement(Transform3::Identity());
  std::shared_ptr<ActsTests::DetectorElementStub> other =
      makeDummyElement(Transform3::Identity());
  // Reverse ID order: the selection order determines the numerical indices.
  std::shared_ptr<AlignableStructure> group =
      std::make_shared<AlignableStructure>(GeometryIdentifier().withVolume(2),
                                           Transform3::Identity());
  std::shared_ptr<AlignableStructure> otherGroup =
      std::make_shared<AlignableStructure>(GeometryIdentifier().withVolume(1),
                                           Transform3::Identity());

  TargetFixture() {
    group->alignmentMask() = AlignmentMask::All;
    group->addSurface(first->surface().getSharedPtr());
    group->addSurface(second->surface().getSharedPtr());
    otherGroup->alignmentMask() = AlignmentMask::Center2;
    otherGroup->addSurface(other->surface().getSharedPtr());
  }
};
}  // namespace

BOOST_FIXTURE_TEST_CASE(ExplicitStructureSelectionValidation, TargetFixture) {
  auto metadata = validateSelection({group}, {});
  BOOST_REQUIRE(metadata.ok());
  BOOST_CHECK(!*metadata);
  auto parent = std::make_shared<AlignableStructure>(GeometryIdentifier());
  parent->addChild(group);
  auto selected = validateSelection({parent}, {group});
  BOOST_REQUIRE(selected.ok());
  BOOST_CHECK(*selected);
  parent->alignmentMask() = AlignmentMask::All;
  checkTargetError(validateSelection({parent}, {parent}),
                   AlignmentError::InvalidAlignmentSelection, "non-leaf");
  checkTargetError(validateSelection({parent}, {parent, group}),
                   AlignmentError::InvalidAlignmentSelection, "non-leaf");
  checkTargetError(validateSelection({group}, {nullptr}),
                   AlignmentError::InvalidAlignmentSelection, "null");
  checkTargetError(validateSelection({group}, {otherGroup}),
                   AlignmentError::InvalidAlignmentSelection, "unregistered");
  checkTargetError(validateSelection({group}, {group, group}),
                   AlignmentError::InvalidAlignmentSelection, "Duplicate");
  parent->alignmentMask() = AlignmentMask::None;
  group->alignmentMask() = AlignmentMask::None;
  group->constraints()[eAlignmentCenter0] = 0.;
  group->setTransform([](const GeometryContext&) -> Transform3 {
    throw std::runtime_error("validation must not query providers");
  });
  auto inactive = validateSelection({parent}, {parent, group});
  BOOST_REQUIRE(inactive.ok());
  BOOST_CHECK(!*inactive);
  checkTargetError(validateSelection({parent}, {group, group}),
                   AlignmentError::InvalidAlignmentSelection, "Duplicate");
}

BOOST_FIXTURE_TEST_CASE(SelectionSnapshotsFramesMasksAndMembership,
                        TargetFixture) {
  const Transform3 initial =
      Translation3(5_mm, -20_mm, 8_mm) * AngleAxis3(0.4, Vector3::UnitX());
  GeometryContext context{initial};
  std::size_t calls = 0;
  group->setTransform([&](const GeometryContext& ctx) {
    ++calls;
    return ctx.get<Transform3>();
  });
  auto prepared = prepareTargets(context, {group}, {group});
  BOOST_REQUIRE(prepared.ok());
  BOOST_CHECK_EQUAL(calls, 1u);
  group->setTransform(Transform3::Identity());
  group->alignmentMask() = AlignmentMask::None;
  group->addSurface(other->surface().getSharedPtr());
  BOOST_CHECK(prepared->targets[0].mask == AlignmentMask::All);
  BOOST_CHECK_EQUAL(prepared->targets[0].members.size(), 2u);
  CHECK_CLOSE_ABS(prepared->targets[0].frame->matrix(), initial.matrix(),
                  1e-14);
  CHECK_CLOSE_ABS(prepared->targets[0].members[1].transform.matrix(),
                  second->localToGlobalTransform(gctx).matrix(), 1e-14);
  std::weak_ptr<const AlignableStructure> lifetime = group;
  group.reset();
  BOOST_CHECK(!lifetime.expired());
}

BOOST_FIXTURE_TEST_CASE(InactiveMetadataAndSelectedLeaves, TargetFixture) {
  auto parent = std::make_shared<AlignableStructure>(GeometryIdentifier());
  parent->addChild(group);
  auto prepared = prepareTargets(gctx, {parent}, {group});
  BOOST_REQUIRE(prepared.ok());
  BOOST_CHECK_EQUAL(prepared->targets.size(), 1u);
  parent->alignmentMask() = AlignmentMask::All;
  checkTargetError(prepareTargets(gctx, {parent}, {parent}),
                   AlignmentError::InvalidAlignmentSelection, "non-leaf");
  checkTargetError(prepareTargets(gctx, {parent}, {parent, group}),
                   AlignmentError::InvalidAlignmentSelection, "non-leaf");
  parent->alignmentMask() = AlignmentMask::None;
  group->alignmentMask() = AlignmentMask::None;
  group->constraints()[eAlignmentCenter0] = 0.;
  group->setTransform([](const GeometryContext&) -> Transform3 {
    throw std::runtime_error("inactive provider must not be queried");
  });
  auto inactive =
      prepareTargets(gctx, {parent}, {parent, group}, {first.get()});
  BOOST_REQUIRE(inactive.ok());
  BOOST_REQUIRE_EQUAL(inactive->targets.size(), 1u);
  BOOST_CHECK(inactive->targets.front().kind ==
              AlignDetail::AlignmentTarget::Kind::Module);
  auto empty = prepareTargets(gctx, {parent}, {group}, {}, {});
  BOOST_REQUIRE(empty.ok());
  BOOST_CHECK(empty->targets.empty());
}

BOOST_FIXTURE_TEST_CASE(InvalidTargetSelections, TargetFixture) {
  checkTargetError(prepareTargets(gctx, {group}, {nullptr}),
                   AlignmentError::InvalidAlignmentSelection, "null");
  checkTargetError(prepareTargets(gctx, {group}, {otherGroup}),
                   AlignmentError::InvalidAlignmentSelection, "unregistered");
  checkTargetError(prepareTargets(gctx, {group}, {group, group}),
                   AlignmentError::InvalidAlignmentSelection, "Duplicate");
  checkTargetError(prepareTargets(gctx, {group}, {group}, {first.get()}),
                   AlignmentError::InvalidAlignmentSelection, "overlapping");
  checkTargetError(prepareTargets(gctx, {}, {}, {nullptr}),
                   AlignmentError::InvalidAlignmentSelection, "Null");
  checkTargetError(prepareTargets(gctx, {}, {}, {module.get(), module.get()}),
                   AlignmentError::InvalidAlignmentSelection, "Duplicate");
  group->alignmentMask() = AlignmentMask::None;
  checkTargetError(prepareTargets(gctx, {group}, {group, group}),
                   AlignmentError::InvalidAlignmentSelection, "Duplicate");
  auto empty = std::make_shared<AlignableStructure>(GeometryIdentifier());
  empty->alignmentMask() = AlignmentMask::All;
  checkTargetError(prepareTargets(gctx, {empty}, {empty}),
                   AlignmentError::InvalidAlignmentSelection, "Empty active");
  empty->addSurface(Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(1_mm, 1_mm)));
  checkTargetError(prepareTargets(gctx, {empty}, {empty}),
                   AlignmentError::InvalidAlignmentSelection, "placement");
}

BOOST_FIXTURE_TEST_CASE(InvalidGraphsAreRejectedBeforeTraversal,
                        TargetFixture) {
  checkTargetError(validateSelection({nullptr}, {}),
                   AlignmentError::HierarchyValidationFailure, "Null");
  checkTargetError(validateSelection({group, group}, {}),
                   AlignmentError::HierarchyValidationFailure, "repeated");
  auto parent = std::make_shared<AlignableStructure>(GeometryIdentifier());
  parent->addChild(group);
  parent->addChild(group);
  checkTargetError(validateSelection({parent}, {}),
                   AlignmentError::HierarchyValidationFailure, "repeated");
  auto nullChild = std::make_shared<AlignableStructure>(GeometryIdentifier());
  nullChild->addChild(nullptr);
  checkTargetError(validateSelection({nullChild}, {}),
                   AlignmentError::HierarchyValidationFailure, "Null");
  auto nullSurface = std::make_shared<AlignableStructure>(GeometryIdentifier());
  nullSurface->addSurface(nullptr);
  checkTargetError(validateSelection({nullSurface}, {}),
                   AlignmentError::HierarchyValidationFailure, "Null member");
  otherGroup->addSurface(first->surface().getSharedPtr());
  checkTargetError(validateSelection({group, otherGroup}, {group, otherGroup}),
                   AlignmentError::HierarchyValidationFailure, "placement");
  group->addSurface(first->surface().getSharedPtr());
  checkTargetError(validateSelection({group}, {group}),
                   AlignmentError::HierarchyValidationFailure, "placement");
  AlignableStructure cycleA{GeometryIdentifier()}, cycleB{GeometryIdentifier()};
  auto a = std::shared_ptr<AlignableStructure>(&cycleA, [](auto*) {});
  auto b = std::shared_ptr<AlignableStructure>(&cycleB, [](auto*) {});
  cycleA.addChild(b);
  cycleB.addChild(a);
  checkTargetError(validateSelection({a}, {}),
                   AlignmentError::HierarchyValidationFailure, "Cycle");
}

BOOST_FIXTURE_TEST_CASE(ConstraintsAndUpdateCapabilities, TargetFixture) {
  for (double variance : {0., 1., -1.}) {
    group->constraints()[eAlignmentCenter0] = variance;
    checkTargetError(prepareTargets(gctx, {group}, {group}),
                     AlignmentError::UnsupportedAlignmentConstraints,
                     "Constraints");
  }
  group->constraints().clear();
  checkTargetError(prepareTargets(gctx, {group}, {group}, {}, {}),
                   AlignmentError::InvalidAlignmentSelection, "update route");
  std::size_t calls = 0;
  checkTargetError(prepareTargets(gctx, {group}, {group}, {},
                                  [&](const auto& placement) {
                                    ++calls;
                                    return &placement != second.get();
                                  }),
                   AlignmentError::InvalidAlignmentSelection, "update route");
  BOOST_CHECK_EQUAL(calls, 2u);
  checkTargetError(
      prepareTargets(gctx, {group}, {group}, {},
                     [](const auto&) -> bool {
                       throw std::runtime_error("capability failed");
                     }),
      AlignmentError::InvalidAlignmentSelection, "capability failed");
}

BOOST_FIXTURE_TEST_CASE(InvalidFramesAndProviderExceptions, TargetFixture) {
  std::vector<Transform3> invalid(6, Transform3::Identity());
  invalid[0].matrix()(0, 3) = std::numeric_limits<double>::quiet_NaN();
  invalid[1].matrix()(0, 0) = std::numeric_limits<double>::infinity();
  invalid[2].linear()(0, 0) = 2.;
  invalid[3].linear()(0, 1) = 0.1;
  invalid[4].linear()(0, 0) = -1.;
  invalid[5].matrix()(3, 0) = 0.1;
  for (const auto& transform : invalid) {
    group->setTransform(transform);
    checkTargetError(prepareTargets(gctx, {group}, {group}),
                     AlignmentError::AlignmentFrameResolutionFailure,
                     "invalid initial structure frame");
  }
  auto noFrame = std::make_shared<AlignableStructure>(GeometryIdentifier());
  noFrame->alignmentMask() = AlignmentMask::All;
  noFrame->addSurface(first->surface().getSharedPtr());
  checkTargetError(prepareTargets(gctx, {noFrame}, {noFrame}),
                   AlignmentError::AlignmentFrameResolutionFailure,
                   "no frame provider");
  group->setTransform([](const GeometryContext&) -> Transform3 {
    throw std::runtime_error("conditions IOV missing");
  });
  auto failure = prepareTargets(gctx, {group}, {group});
  checkTargetError(failure, AlignmentError::AlignmentFrameResolutionFailure,
                   "conditions IOV missing");
  BOOST_CHECK_NE(failure.error().message.find("fit geometry context"),
                 std::string::npos);
  group->setTransform(
      [](const GeometryContext&) -> Transform3 { throw std::bad_alloc(); });
  BOOST_CHECK_THROW(static_cast<void>(prepareTargets(gctx, {group}, {group})),
                    std::bad_alloc);
  group->setTransform([](const GeometryContext&) -> Transform3 { throw 42; });
  BOOST_CHECK_THROW(static_cast<void>(prepareTargets(gctx, {group}, {group})),
                    int);
}

namespace {
class ContextualTargetElement : public ActsTests::DetectorElementStub {
 public:
  ContextualTargetElement()
      : DetectorElementStub(
            Transform3::Identity(),
            std::make_shared<const RectangleBounds>(10_mm, 10_mm), 1_um) {}

  const Transform3& localToGlobalTransform(
      const GeometryContext& ctx) const override {
    ++calls;
    return ctx.get<Transform3>();
  }
  mutable std::size_t calls = 0;
};

class WrongBacklinkElement : public ActsTests::DetectorElementStub {
 public:
  explicit WrongBacklinkElement(const Surface& other)
      : DetectorElementStub(
            Transform3::Identity(),
            std::make_shared<const RectangleBounds>(10_mm, 10_mm), 1_um),
        m_other(other) {}
  using DetectorElementStub::surface;
  const Surface& surface() const override { return m_other; }

 private:
  const Surface& m_other;
};
}  // namespace

BOOST_FIXTURE_TEST_CASE(MemberContextAndBacklinkValidation, TargetFixture) {
  ContextualTargetElement contextual;
  auto selected = std::make_shared<AlignableStructure>(GeometryIdentifier(),
                                                       Transform3::Identity());
  selected->alignmentMask() = AlignmentMask::All;
  selected->addSurface(contextual.surface().getSharedPtr());
  const Transform3 initial =
      Translation3(10_mm, -40_mm, 5_mm) * AngleAxis3(0.2, Vector3::UnitZ());
  GeometryContext context{initial};
  auto prepared = prepareTargets(context, {selected}, {selected});
  BOOST_REQUIRE(prepared.ok());
  BOOST_CHECK_EQUAL(contextual.calls, 1u);
  CHECK_CLOSE_ABS(prepared->targets[0].members[0].transform.matrix(),
                  initial.matrix(), 1e-14);
  checkTargetError(prepareTargets(gctx, {selected}, {selected}),
                   AlignmentError::AlignmentFrameResolutionFailure,
                   "initial member 0 frame in fit geometry context");
  Transform3 invalid = initial;
  invalid.matrix()(3, 3) = 0.;
  GeometryContext invalidContext{invalid};
  checkTargetError(prepareTargets(invalidContext, {selected}, {selected}),
                   AlignmentError::AlignmentFrameResolutionFailure,
                   "invalid initial member");
  WrongBacklinkElement wrong(first->surface());
  checkTargetError(prepareTargets(gctx, {}, {}, {&wrong}),
                   AlignmentError::InvalidAlignmentSelection, "backlink");
  auto wrongGroup = std::make_shared<AlignableStructure>(
      GeometryIdentifier(), Transform3::Identity());
  wrongGroup->alignmentMask() = AlignmentMask::Center0;
  wrongGroup->addSurface(wrong.surface().getSharedPtr());
  checkTargetError(prepareTargets(gctx, {wrongGroup}, {wrongGroup}),
                   AlignmentError::InvalidAlignmentSelection, "backlink");
}

BOOST_FIXTURE_TEST_CASE(PreparedTargetOrdering, TargetFixture) {
  auto prepared = prepareTargets(gctx, {otherGroup, group}, {group, otherGroup},
                                 {module.get()});
  BOOST_REQUIRE(prepared.ok());
  BOOST_REQUIRE_EQUAL(prepared->targets.size(), 3u);
  const auto& standalone = prepared->targets[0];
  BOOST_CHECK(standalone.kind == AlignDetail::AlignmentTarget::Kind::Module);
  BOOST_CHECK(!standalone.structure);
  BOOST_CHECK(!standalone.frame);
  BOOST_REQUIRE_EQUAL(standalone.members.size(), 1u);
  BOOST_CHECK(standalone.members[0].placement == module.get());
  CHECK_CLOSE_ABS(standalone.members[0].transform.matrix(),
                  module->localToGlobalTransform(gctx).matrix(), 1e-14);
  BOOST_CHECK(prepared->targets[1].kind ==
              AlignDetail::AlignmentTarget::Kind::Structure);
  BOOST_CHECK(prepared->targets[1].structure == group);
  BOOST_CHECK(prepared->targets[2].structure == otherGroup);
  BOOST_CHECK(prepared->targets[2].mask == AlignmentMask::Center2);
}

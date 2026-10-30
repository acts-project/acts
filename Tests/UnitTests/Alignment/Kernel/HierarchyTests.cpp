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
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsAlignment/Geometry/AlignableStructure.hpp"
#include "ActsAlignment/Geometry/AlignmentHierarchy.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"

#include <memory>
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

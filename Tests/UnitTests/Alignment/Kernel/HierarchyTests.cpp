// This file is part of the ACTS project.
//
// Copyright (C) 2026 CERN for the benefit of the ACTS project
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
  Acts::Result<HierarchyValidationTestTrack> fit(it_t /*begin*/, it_t /*end*/,
                                                 const start_parameters_t&,
                                                 const fit_options_t&,
                                                 track_container_t&) const {
    return Acts::Result<HierarchyValidationTestTrack>::success(
        HierarchyValidationTestTrack{});
  }
};

}  // namespace

BOOST_AUTO_TEST_CASE(HierarchyValidation) {
  const auto geoCtx = Acts::GeometryContext::dangerouslyDefaultConstruct();

  // 1. Create 3 detector elements
  auto el1 = makeDummyElement(Acts::Translation3(0, 0, 10_mm) *
                              Acts::Transform3::Identity());
  auto el2 = makeDummyElement(Acts::Translation3(0, 0, 20_mm) *
                              Acts::Transform3::Identity());
  auto el3 = makeDummyElement(Acts::Translation3(0, 0, 30_mm) *
                              Acts::Transform3::Identity());

  std::vector<Acts::SurfacePlacementBase*> elements = {el1.get(), el2.get(),
                                                       el3.get()};

  // 2. Setup structures
  auto struct1 = std::make_shared<ActsAlignment::AlignableStructure>(
      Acts::GeometryIdentifier().withVolume(1));
  auto struct2 = std::make_shared<ActsAlignment::AlignableStructure>(
      Acts::GeometryIdentifier().withVolume(2));

  // Dummy alignment options
  DummyFitter fitter;
  ActsAlignment::Alignment alignEngine(std::move(fitter));

  ActsAlignment::AlignedTransformUpdater voidUpdater =
      [](Acts::SurfacePlacementBase*, const Acts::GeometryContext&,
         const Acts::Transform3&) { return true; };

  DummyFitOptions kfOptions;
  kfOptions.geoContext = geoCtx;

  // --- Test Case 1: Mixed mode ---
  // Only el1 and el2 are in structures. el3 floats as a standalone module.
  // Mixed mode is allowed: the hierarchy check must not flag this.
  struct1->addSurface(el1->surface().getSharedPtr());
  struct2->addSurface(el2->surface().getSharedPtr());

  ActsAlignment::AlignmentOptions<DummyFitOptions> options(
      kfOptions, voidUpdater, elements, 0.5, {5, 0.01}, 5, {},
      {struct1.get(), struct2.get()});

  // Create dummy collections for the align call
  std::vector<std::vector<Acts::SourceLink>> trajs;
  std::vector<Acts::BoundTrackParameters> params;

  auto res1 = alignEngine.align(trajs, params, options);
  if (!res1.ok()) {
    BOOST_CHECK_NE(res1.error(),
                   ActsAlignment::AlignmentError::HierarchyValidationFailure);
  }

  // --- Test Case 2: Overlap detection ---
  // Add el1 to struct2 as well, so el1 is now in both struct1 and struct2.
  struct2->addSurface(el1->surface().getSharedPtr());

  ActsAlignment::AlignmentOptions<DummyFitOptions> options2(
      kfOptions, voidUpdater, elements, 0.5, {5, 0.01}, 5, {},
      {struct1.get(), struct2.get()});
  auto res2 = alignEngine.align(trajs, params, options2);
  BOOST_CHECK(!res2.ok());
  BOOST_CHECK_EQUAL(res2.error(),
                    ActsAlignment::AlignmentError::HierarchyValidationFailure);

  // --- Test Case 3: Correct assignment ---
  // Clear and re-assign correctly
  auto structA = std::make_shared<ActsAlignment::AlignableStructure>(
      Acts::GeometryIdentifier().withVolume(10));
  auto structB = std::make_shared<ActsAlignment::AlignableStructure>(
      Acts::GeometryIdentifier().withVolume(20));

  structA->addSurface(el1->surface().getSharedPtr());
  structA->addSurface(el2->surface().getSharedPtr());
  structB->addSurface(el3->surface().getSharedPtr());

  ActsAlignment::AlignmentOptions<DummyFitOptions> options3(
      kfOptions, voidUpdater, elements, 0.5, {5, 0.01}, 5, {},
      {structA.get(), structB.get()});

  auto res3 = alignEngine.align(trajs, params, options3);
  // It will likely fail with something else because we have 0 tracks,
  // but it should PASS the hierarchy check.
  if (!res3.ok()) {
    BOOST_CHECK_NE(res3.error(),
                   ActsAlignment::AlignmentError::HierarchyValidationFailure);
  }
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
    ActsAlignment::AlignmentHierarchy hierarchy(
        {structA.get(), structB.get()});
    BOOST_CHECK(hierarchy.validate().ok());
    BOOST_CHECK_EQUAL(hierarchy.structureFor(el1.get()), structA.get());
    BOOST_CHECK_EQUAL(hierarchy.structureFor(el2.get()), structA.get());
    // el3 is a standalone floating module
    BOOST_CHECK_EQUAL(hierarchy.structureFor(el3.get()), nullptr);
  }

  // --- Overlap detected ---
  structB->addSurface(el1->surface().getSharedPtr());
  {
    ActsAlignment::AlignmentHierarchy hierarchy(
        {structA.get(), structB.get()});
    const auto result = hierarchy.validate();
    BOOST_CHECK(!result.ok());
    BOOST_CHECK_EQUAL(result.overlapping.size(), 1u);
    BOOST_CHECK_EQUAL(result.overlapping.front(), el1.get());
  }

  // --- Empty hierarchy ---
  {
    ActsAlignment::AlignmentHierarchy hierarchy({});
    BOOST_CHECK(hierarchy.validate().ok());
    BOOST_CHECK_EQUAL(hierarchy.structureFor(el1.get()), nullptr);
  }
}

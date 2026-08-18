// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Geant4/AlgebraConverters.hpp"
#include "ActsExamples/Geant4/Geant4ConstructionOptions.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeG4DetectorConstruction.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4RotationMatrix.hh>
#include <G4VPhysicalVolume.hh>

using namespace ActsExamples;

namespace {

/// A sensitive Geant4 volume in global coordinates
struct G4Sensitive {
  Acts::Vector3 center;
  Acts::Vector2 halfLengths;
};

/// Collect all silicon volumes of the tree, composing transforms the way
/// `SensitiveSurfaceMapper` does
void collectSensitives(const G4VPhysicalVolume& physicalVolume,
                       const Acts::Transform3& motherTransform,
                       std::vector<G4Sensitive>& sensitives) {
  const G4LogicalVolume& logicalVolume = *physicalVolume.GetLogicalVolume();

  Acts::Transform3 localToGlobal =
      motherTransform * Acts::Translation3(Geant4::convertPosition(
                            physicalVolume.GetTranslation()));
  if (const G4RotationMatrix* g4Rotation = physicalVolume.GetRotation();
      g4Rotation != nullptr) {
    // Geant4 stores the rotation of the mother relative to the daughter
    // frame, hence the transpose
    Acts::RotationMatrix3 rotation;
    rotation << g4Rotation->xx(), g4Rotation->yx(), g4Rotation->zx(),
        g4Rotation->xy(), g4Rotation->yy(), g4Rotation->zy(), g4Rotation->xz(),
        g4Rotation->yz(), g4Rotation->zz();
    localToGlobal.rotate(rotation);
  }

  if (logicalVolume.GetMaterial()->GetName() == "Silicon") {
    const auto& box = dynamic_cast<const G4Box&>(*logicalVolume.GetSolid());
    constexpr double convertLength = CLHEP::mm / Acts::UnitConstants::mm;
    sensitives.push_back({localToGlobal.translation(),
                          {box.GetXHalfLength() * convertLength,
                           box.GetYHalfLength() * convertLength}});
  }

  for (std::size_t i = 0; i < logicalVolume.GetNoDaughters(); ++i) {
    collectSensitives(*logicalVolume.GetDaughter(i), localToGlobal, sensitives);
  }
}

/// Check the whole tree for overlaps
bool hasOverlaps(G4VPhysicalVolume& physicalVolume) {
  constexpr G4int resolution = 1000;
  constexpr G4double tolerance = 0.;
  constexpr G4bool verbose = false;
  bool overlaps = physicalVolume.CheckOverlaps(resolution, tolerance, verbose);

  G4LogicalVolume& logicalVolume = *physicalVolume.GetLogicalVolume();
  for (std::size_t i = 0; i < logicalVolume.GetNoDaughters(); ++i) {
    overlaps |= hasOverlaps(*logicalVolume.GetDaughter(i));
  }
  return overlaps;
}

TelescopeDetector::Config makeConfig(int binValue) {
  TelescopeDetector::Config cfg;
  cfg.positions = {30, 60, 90, 120, 150, 180};
  cfg.stereos = {0, 0, 0, 0, 0, 0};
  cfg.offsets = {10, -20};
  cfg.bounds = {25, 100};
  cfg.binValue = binValue;
  return cfg;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(TelescopeG4DetectorConstructionTests)

BOOST_AUTO_TEST_CASE(ConstructSmoke) {
  TelescopeDetector detector{TelescopeDetector::Config{}};

  auto construction =
      detector.buildGeant4DetectorConstruction(Geant4ConstructionOptions{});
  BOOST_REQUIRE(construction != nullptr);

  G4VPhysicalVolume* world = construction->Construct();
  BOOST_REQUIRE(world != nullptr);
  // The world is built once and cached
  BOOST_CHECK_EQUAL(construction->Construct(), world);

  BOOST_CHECK(!hasOverlaps(*world));
}

// The Geant4 sensors have to match the tracking geometry surfaces on all axes
BOOST_DATA_TEST_CASE(MatchesTrackingGeometry,
                     boost::unit_test::data::make({0, 1, 2}), binValue) {
  const auto cfg = makeConfig(binValue);
  TelescopeDetector detector{cfg};

  auto construction =
      detector.buildGeant4DetectorConstruction(Geant4ConstructionOptions{});
  G4VPhysicalVolume* world = construction->Construct();
  BOOST_REQUIRE(world != nullptr);

  BOOST_CHECK(!hasOverlaps(*world));

  std::vector<G4Sensitive> sensitives;
  collectSensitives(*world, Acts::Transform3::Identity(), sensitives);
  BOOST_REQUIRE_EQUAL(sensitives.size(), cfg.positions.size());

  const auto& gctx = detector.nominalGeometryContext();
  std::vector<const Acts::Surface*> surfaces;
  detector.trackingGeometry()->visitSurfaces(
      [&](const Acts::Surface* surface) { surfaces.push_back(surface); }, true);
  BOOST_REQUIRE_EQUAL(surfaces.size(), cfg.positions.size());

  for (const auto* surface : surfaces) {
    const Acts::Vector3 center = surface->center(gctx);

    auto it = std::ranges::find_if(sensitives, [&](const G4Sensitive& g4) {
      return (g4.center - center).norm() < 1e-6;
    });
    BOOST_REQUIRE_MESSAGE(
        it != sensitives.end(),
        "no Geant4 volume at surface center " << center.transpose());

    const auto& bounds =
        dynamic_cast<const Acts::RectangleBounds&>(surface->bounds());
    BOOST_CHECK_CLOSE(it->halfLengths[0], bounds.halfLengthX(), 1e-6);
    BOOST_CHECK_CLOSE(it->halfLengths[1], bounds.halfLengthY(), 1e-6);
  }
}

// Unsorted positions and a stack that is not centered on the origin
BOOST_AUTO_TEST_CASE(UnsortedAndShiftedPositions) {
  TelescopeDetector::Config cfg;
  cfg.positions = {-100, 200, 50};
  cfg.stereos = {0, 0, 0};
  cfg.offsets = {40, 40};

  TelescopeDetector detector{cfg};
  auto construction =
      detector.buildGeant4DetectorConstruction(Geant4ConstructionOptions{});
  G4VPhysicalVolume* world = construction->Construct();
  BOOST_REQUIRE(world != nullptr);

  BOOST_CHECK(!hasOverlaps(*world));

  std::vector<G4Sensitive> sensitives;
  collectSensitives(*world, Acts::Transform3::Identity(), sensitives);
  BOOST_CHECK_EQUAL(sensitives.size(), cfg.positions.size());
}

BOOST_AUTO_TEST_CASE(EmptyPositionsThrows) {
  TelescopeDetector::Config cfg;
  cfg.positions = {};
  cfg.stereos = {};

  BOOST_CHECK_THROW(TelescopeDetector{cfg}, std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

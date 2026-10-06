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
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"

#include <algorithm>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4Tubs.hh>
#include <G4VPhysicalVolume.hh>

using namespace Acts::UnitLiterals;
using namespace ActsExamples::Geant4;

namespace {

auto gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();

/// Default translation tolerance of `Acts::detail::TransformComparator`,
/// which the mapper uses to compare surface centers and Geant4 positions
constexpr double kTolerance = 0.1_um;

/// All surfaces are candidates for every position. With `cacheByRegion` they
/// form a single region, so the mapper builds its lookup once; without it the
/// mapper queries the candidates for every Geant4 volume.
struct ListCandidates : public SensitiveCandidatesBase {
  std::vector<const Acts::Surface*> surfaces;
  bool cacheByRegion = false;

  std::vector<const Acts::Surface*> queryPosition(
      const Acts::GeometryContext& /*gctx*/,
      const Acts::Vector3& /*position*/) const override {
    return surfaces;
  }

  std::vector<const Acts::Surface*> queryAll() const override {
    return surfaces;
  }

  std::optional<const Acts::GeometryObject*> queryRegion(
      const Acts::GeometryContext& /*gctx*/,
      const Acts::Vector3& /*position*/) const override {
    if (!cacheByRegion) {
      return std::nullopt;
    }
    return nullptr;
  }
};

std::shared_ptr<Acts::Surface> makePlane(const Acts::Vector3& center,
                                         double halfLength = 0.5_mm) {
  return Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3(Acts::Translation3(center)),
      std::make_shared<Acts::RectangleBounds>(halfLength, halfLength));
}

/// A Geant4 world with daughters placed without rotation
struct World {
  G4Material* material =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
  G4LogicalVolume* logical = new G4LogicalVolume(
      new G4Box("World", 10_m, 10_m, 10_m), material, "World");
  G4VPhysicalVolume* physical = new G4PVPlacement(
      nullptr, G4ThreeVector(), logical, "World", nullptr, false, 0);

  const G4VPhysicalVolume* place(G4VSolid* solid, const std::string& name,
                                 const Acts::Vector3& position) {
    auto* volume = new G4LogicalVolume(solid, material, name);
    return new G4PVPlacement(
        nullptr, G4ThreeVector(position.x(), position.y(), position.z()),
        volume, name, logical, false, 0);
  }

  const G4VPhysicalVolume* placeBox(const std::string& name,
                                    const Acts::Vector3& position,
                                    double halfLength = 0.5_mm) {
    return place(new G4Box(name, halfLength, halfLength, 10_um), name,
                 position);
  }
};

struct MappingResult {
  SensitiveSurfaceMapper::State state;
  std::shared_ptr<SensitiveSurfaceMapper> mapper;

  /// The surface a Geant4 volume is mapped to, nullptr if none
  const Acts::Surface* surfaceOf(const G4VPhysicalVolume* volume) const {
    const auto it = state.g4VolumeToSurfaces.find(volume);
    if (it == state.g4VolumeToSurfaces.end()) {
      return nullptr;
    }
    BOOST_REQUIRE_EQUAL(it->second.size(), 1u);
    return it->second.begin()->second;
  }
};

MappingResult runMapper(World& world,
                        const std::vector<const Acts::Surface*>& surfaces,
                        bool cacheByRegion) {
  auto candidates = std::make_shared<ListCandidates>();
  candidates->surfaces = surfaces;
  candidates->cacheByRegion = cacheByRegion;

  SensitiveSurfaceMapper::Config cfg;
  cfg.volumeMappings = {"Sensor"};
  cfg.candidateSurfaces = candidates;

  MappingResult result;
  result.mapper = std::make_shared<SensitiveSurfaceMapper>(
      cfg,
      Acts::getDefaultLogger("SensitiveSurfaceMapper", Acts::Logging::WARNING));
  result.mapper->remapSensitiveNames(result.state, gctx, world.physical,
                                     Acts::Transform3::Identity());
  return result;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(Geant4SensitiveSurfaceMapperSuite)

BOOST_DATA_TEST_CASE(CenterMatchingTolerance,
                     boost::unit_test::data::make({false, true}),
                     cacheByRegion) {
  World world;
  // Centers on and next to multiples of 1 mm
  const std::vector<Acts::Vector3> centers = {
      {0., 0., 0.}, {1., 2., -3.}, {1.5, 2.5, -3.5}, {-7.99997, 4., 12.}};
  std::vector<std::shared_ptr<Acts::Surface>> planes;
  std::vector<const Acts::Surface*> surfaces;
  for (const auto& center : centers) {
    planes.push_back(makePlane(center));
    surfaces.push_back(planes.back().get());
  }

  std::vector<const G4VPhysicalVolume*> inside;
  std::vector<const G4VPhysicalVolume*> outside;
  for (const auto& center : centers) {
    // Within the tolerance in every component. For the second and last
    // center this crosses a multiple of 1 mm, i.e. a cell of the mapper
    inside.push_back(world.placeBox(
        "Sensor", center + Acts::Vector3(-0.5, 0.5, -0.9) * kTolerance));
    outside.push_back(world.placeBox(
        "Sensor", center + Acts::Vector3(0., 0., 2. * kTolerance)));
  }

  const auto result = runMapper(world, surfaces, cacheByRegion);
  for (std::size_t i = 0; i < centers.size(); ++i) {
    BOOST_CHECK_EQUAL(result.surfaceOf(inside[i]), surfaces[i]);
    BOOST_CHECK_EQUAL(result.surfaceOf(outside[i]), nullptr);
  }
  BOOST_CHECK_EQUAL(result.state.missingVolumes.size(), centers.size());
}

BOOST_DATA_TEST_CASE(SharedCenterConcentricCylinders,
                     boost::unit_test::data::make({false, true}) *
                         boost::unit_test::data::make({false, true}),
                     cacheByRegion, reverseOrder) {
  // Concentric cylinders around the beam line, like a vertex detector, all
  // with the center of their Geant4 tubes
  World world;
  const std::vector<double> radii = {5_mm, 12_mm, 25_mm};
  std::vector<std::shared_ptr<Acts::Surface>> cylinders;
  std::vector<const G4VPhysicalVolume*> tubes;
  for (const double r : radii) {
    cylinders.push_back(Acts::Surface::makeShared<Acts::CylinderSurface>(
        Acts::Transform3::Identity(),
        std::make_shared<Acts::CylinderBounds>(r + 10_um, 250_mm)));
    tubes.push_back(world.place(
        new G4Tubs("Tube", r, r + 20_um, 250_mm, 0., 2. * std::numbers::pi),
        "Sensor", Acts::Vector3::Zero()));
  }
  std::vector<const Acts::Surface*> surfaces;
  for (const auto& cylinder : cylinders) {
    surfaces.push_back(cylinder.get());
  }
  if (reverseOrder) {
    std::ranges::reverse(surfaces);
  }

  const auto result = runMapper(world, surfaces, cacheByRegion);
  for (std::size_t i = 0; i < radii.size(); ++i) {
    BOOST_CHECK_EQUAL(result.surfaceOf(tubes[i]), cylinders[i].get());
  }
}

BOOST_DATA_TEST_CASE(IdenticalSurfacesKeepCandidateOrder,
                     boost::unit_test::data::make({false, true}),
                     cacheByRegion) {
  World world;
  const Acts::Vector3 center{10., -20., 30.};
  const auto first = makePlane(center);
  const auto second = makePlane(center);
  const auto* box = world.placeBox("Sensor", center);

  const auto result =
      runMapper(world, {first.get(), second.get()}, cacheByRegion);
  BOOST_CHECK_EQUAL(result.surfaceOf(box), first.get());
}

BOOST_DATA_TEST_CASE(AnnulusCentroidAndCenterMatchFollowCandidateOrder,
                     boost::unit_test::data::make({false, true}) *
                         boost::unit_test::data::make({false, true}),
                     cacheByRegion, annulusFirst) {
  // An annulus whose bounds centroid lies inside a large Geant4 box, and a
  // plane centered on the same box: whichever comes first in the candidate
  // order is the match
  World world;
  const Acts::Vector3 discCenter{0., 0., 100.};
  const auto annulus = Acts::Surface::makeShared<Acts::DiscSurface>(
      Acts::Transform3(Acts::Translation3(discCenter)),
      std::make_shared<Acts::AnnulusBounds>(380_mm, 500_mm, -0.07, 0.07,
                                            Acts::Vector2{2_mm, -1_mm}));
  const Acts::Vector3 boxCenter = discCenter + Acts::Vector3(440., 0., 0.);
  const auto plane = makePlane(boxCenter);
  const auto* box =
      world.place(new G4Box("BigBox", 1_m, 1_m, 1_m), "Sensor", boxCenter);

  std::vector<const Acts::Surface*> surfaces = {plane.get(), annulus.get()};
  if (annulusFirst) {
    std::ranges::reverse(surfaces);
  }
  const auto result = runMapper(world, surfaces, cacheByRegion);
  BOOST_CHECK_EQUAL(result.surfaceOf(box), surfaces.front());
}

BOOST_AUTO_TEST_CASE(CheckMappingRequiresOnlySurfacesFlaggedSensitive) {
  World world;
  const auto sensor = makePlane({0., 0., 0.});
  sensor->assignIsSensitive(true);
  // E.g. a passive surface with a sensitive geometry id, without Geant4
  // counterpart
  const auto passive = makePlane({0., 0., 50.});
  passive->assignIsSensitive(false);
  world.placeBox("Sensor", {0., 0., 0.});

  const auto result = runMapper(world, {sensor.get(), passive.get()}, true);
  BOOST_CHECK(result.mapper->checkMapping(result.state, gctx));

  // A sensor without Geant4 counterpart is still reported
  passive->assignIsSensitive(true);
  BOOST_CHECK(!result.mapper->checkMapping(result.state, gctx));
}

BOOST_AUTO_TEST_SUITE_END()

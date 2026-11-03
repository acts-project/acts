// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Material/TrackingGeometryMaterial.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"

using namespace Acts;

namespace {
std::shared_ptr<CylinderSurface> keyedSurface(std::uint64_t id,
                                              std::string key) {
  auto surface =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20., 100.);
  surface->assignGeometryId(GeometryIdentifier().withSensitive(id));
  surface->assignSurfaceMaterial(std::make_shared<ProtoSurfaceMaterial>(
      BinUtility{}, MappingType::Default, std::move(key)));
  return surface;
}

std::shared_ptr<const ISurfaceMaterial> material(double thickness) {
  return std::make_shared<HomogeneousSurfaceMaterial>(
      MaterialSlab(Material::fromMolarDensity(1., 2., 3., 4., 5.), thickness));
}
}  // namespace

BOOST_AUTO_TEST_SUITE(TrackingGeometryMaterialTests)

BOOST_AUTO_TEST_CASE(ExtractPreservesMaterialAndDeferredBinning) {
  auto world = std::make_shared<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(0., 100., 200.), "world");
  world->assignGeometryId(GeometryIdentifier().withVolume(1));
  auto volumeMaterial =
      std::make_shared<HomogeneousVolumeMaterial>(Material::Vacuum());
  world->assignVolumeMaterial(volumeMaterial);
  auto surface =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30., 100.);
  surface->assignGeometryId(GeometryIdentifier().withSensitive(1));
  auto proto = std::make_shared<ProtoGridSurfaceMaterial>(
      MultiAxisSpec2D(
          {AxisSpec::DeferredEquidistant(8, AxisDirection::AxisRPhi),
           AxisSpec::DeferredEquidistant(3, AxisDirection::AxisZ)}),
      MappingType::Default, "cylinder");
  surface->assignSurfaceMaterial(proto);
  world->addSurface(surface);
  auto bare =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 40., 100.);
  bare->assignGeometryId(GeometryIdentifier().withSensitive(2));
  world->addSurface(bare);
  TrackingGeometry geometry(world, nullptr, {}, getDummyLogger(), false);

  const auto material = TrackingGeometryMaterial::fromGeometry(geometry);
  BOOST_CHECK(material.surfaceMaterials.empty());
  BOOST_REQUIRE_EQUAL(material.keyedSurfaces.size(), 1);
  BOOST_CHECK(material.keyedSurfaces.at("cylinder").material == proto);
  BOOST_CHECK(material.volumeMaterials.at(world->geometryId()) ==
              volumeMaterial);

  const auto all = TrackingGeometryMaterial::fromGeometry(geometry, true);
  BOOST_CHECK(all.keyedSurfaces.at("cylinder").material == proto);
  const auto* placeholder = dynamic_cast<const ProtoGridSurfaceMaterial*>(
      all.surfaceMaterials.at(bare->geometryId()).get());
  BOOST_REQUIRE(placeholder != nullptr);
  BOOST_CHECK(placeholder->binning().isDeferred());
  for (const auto& axis : placeholder->binning().axisSpecs()) {
    BOOST_CHECK_EQUAL(axis.asEquidistant().nBins, 1);
  }
  BOOST_CHECK(bare->surfaceMaterial() == nullptr);
  BOOST_CHECK(surface->surfaceMaterialSharedPtr() == proto);
  BOOST_CHECK(proto->binning().isDeferred());
}

BOOST_AUTO_TEST_CASE(ApplyGeometryResolvesKeysAndVolumesBeforeMutation) {
  auto world = std::make_shared<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(0., 100., 200.), "world");
  world->assignGeometryId(GeometryIdentifier().withVolume(1));
  auto a = keyedSurface(2, "a");
  auto b = keyedSurface(1, "b");
  auto legacy =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30., 100.);
  legacy->assignGeometryId(GeometryIdentifier().withSensitive(3));
  world->addPortal(
      std::make_shared<Portal>(GeometryContext::dangerouslyDefaultConstruct(),
                               Portal::Arguments{.alongNormal = {a, *world}}));
  world->addSurface(b);
  world->addSurface(legacy);
  TrackingGeometry geometry(world, nullptr, {}, getDummyLogger(), false);
  auto originalA = a->surfaceMaterialSharedPtr();
  auto originalB = b->surfaceMaterialSharedPtr();
  auto matA = material(1.);
  auto matB = material(2.);
  auto volumeMaterial =
      std::make_shared<HomogeneousVolumeMaterial>(Material::Vacuum());
  TrackingGeometryMaterial maps;
  maps.keyedSurfaces.emplace("a", KeyedSurfaceMaterial{b->geometryId(), matA});
  maps.keyedSurfaces.emplace("unused",
                             KeyedSurfaceMaterial{a->geometryId(), matB});
  // A matching ID must not mask a missing key.
  maps.surfaceMaterials.emplace(b->geometryId(), matA);
  maps.surfaceMaterials.emplace(legacy->geometryId(), matB);
  maps.volumeMaterials.emplace(world->geometryId(), volumeMaterial);
  BOOST_CHECK_THROW(maps.apply(geometry), std::invalid_argument);
  BOOST_CHECK(a->surfaceMaterialSharedPtr() == originalA);
  BOOST_CHECK(b->surfaceMaterialSharedPtr() == originalB);
  BOOST_CHECK(legacy->surfaceMaterial() == nullptr);
  BOOST_CHECK(world->volumeMaterial() == nullptr);

  maps.keyedSurfaces.emplace("b", KeyedSurfaceMaterial{a->geometryId(), matB});
  maps.apply(geometry);
  BOOST_CHECK(a->surfaceMaterialSharedPtr() == matA);
  BOOST_CHECK(b->surfaceMaterialSharedPtr() == matB);
  BOOST_CHECK(legacy->surfaceMaterialSharedPtr() == matB);
  BOOST_CHECK(world->volumeMaterial() == volumeMaterial.get());
}

BOOST_AUTO_TEST_CASE(KeyedMapsRejectGen1BeforeMutation) {
  auto world = std::make_shared<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(0., 100., 200.), "world");
  TrackingGeometry geometry(world);
  BOOST_REQUIRE(geometry.geometryVersion() ==
                TrackingGeometry::GeometryVersion::Gen1);
  TrackingGeometryMaterial maps;
  maps.volumeMaterials.emplace(
      world->geometryId(),
      std::make_shared<HomogeneousVolumeMaterial>(Material::Vacuum()));
  maps.keyedSurfaces.emplace("unused", KeyedSurfaceMaterial{{}, material(1.)});
  BOOST_CHECK_EXCEPTION(
      maps.apply(geometry), std::invalid_argument, [](const auto& error) {
        return std::string(error.what()).find("keyed material map to Gen1") !=
               std::string::npos;
      });
  BOOST_CHECK(world->volumeMaterial() == nullptr);
  maps.keyedSurfaces.clear();
  BOOST_CHECK_NO_THROW(maps.apply(geometry));
  BOOST_CHECK(world->volumeMaterial() != nullptr);
}

BOOST_AUTO_TEST_CASE(SelectedSurfacesRejectAmbiguityAndDeduplicateVisits) {
  auto a = keyedSurface(1, "a");
  auto b = keyedSurface(2, "a");
  auto mat = material(1.);
  TrackingGeometryMaterial maps;
  maps.keyedSurfaces.emplace("a", KeyedSurfaceMaterial{a->geometryId(), mat});
  auto before = a->surfaceMaterialSharedPtr();
  std::vector<Surface*> surfaces{a.get(), b.get()};
  BOOST_CHECK_THROW(maps.apply(surfaces), std::invalid_argument);
  BOOST_CHECK(a->surfaceMaterialSharedPtr() == before);
  surfaces = {a.get(), a.get()};
  BOOST_CHECK_NO_THROW(maps.apply(surfaces));
  BOOST_CHECK(a->surfaceMaterialSharedPtr() == mat);
}

BOOST_AUTO_TEST_SUITE_END()

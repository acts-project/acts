// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/BinnedSurfaceMaterialAccumulator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Material/MaterialMapper.hpp"
#include "Acts/Material/MergedMaterialMarker.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Material/detail/MaterialSurfaceRegistry.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/Diagnostics.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Json/IVolumeMaterialJsonDecorator.hpp"
#include "ActsPlugins/Json/JsonMaterialDecorator.hpp"
#include "ActsPlugins/Json/MaterialMapJsonConverter.hpp"
#include "ActsPlugins/Json/SurfaceJsonConverter.hpp"
#include "ActsPlugins/Json/SurfaceMaterialJsonConverter.hpp"
#include "ActsTests/CommonHelpers/DataDirectory.hpp"

#include <cstdio>
#include <fstream>
#include <memory>

#include <nlohmann/json.hpp>

// These tests intentionally exercise the deprecated material format.
ACTS_PUSH_IGNORE_DEPRECATED()

namespace Acts {
class IVolumeMaterial;
}  // namespace Acts

using namespace Acts;

class DummyDecorator : public IVolumeMaterialJsonDecorator {
 public:
  void decorate([[maybe_unused]] const ISurfaceMaterial& material,
                [[maybe_unused]] nlohmann::json& json) const override {};

  void decorate([[maybe_unused]] const IVolumeMaterial& material,
                [[maybe_unused]] nlohmann::json& json) const override {};
};

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(RoundtripFromFile) {
  // read reference map from file
  std::ifstream refFile(ActsTests::getDataPath("material-map.json"));
  nlohmann::json refJson;
  refFile >> refJson;

  DummyDecorator decorator;
  // convert to the material map and back again
  MaterialMapJsonConverter::Config converterCfg;
  MaterialMapJsonConverter converter(converterCfg, Logging::INFO);
  auto materialMap = converter.jsonToMaterialMaps(refJson);
  nlohmann::json encodedJson =
      converter.materialMapsToJson(materialMap, &decorator);

  // verify identical encoded JSON values
  BOOST_CHECK_EQUAL(refJson, encodedJson);
}

namespace {
std::shared_ptr<CylinderSurface> keyedCylinder(std::uint64_t id,
                                               const std::string& key,
                                               double radius = 20.) {
  auto surface = Surface::makeShared<CylinderSurface>(Transform3::Identity(),
                                                      radius, 100.);
  surface->assignGeometryId(GeometryIdentifier().withSensitive(id));
  surface->assignSurfaceMaterial(std::make_shared<ProtoGridSurfaceMaterial>(
      MultiAxisSpec2D(
          {AxisSpec::DeferredEquidistant(2, AxisDirection::AxisRPhi),
           AxisSpec::DeferredEquidistant(2, AxisDirection::AxisZ)}),
      MappingType::Default, key));
  return surface;
}

std::shared_ptr<const ISurfaceMaterial> materialWithThickness(
    double thickness) {
  return std::make_shared<HomogeneousSurfaceMaterial>(
      MaterialSlab(Material::fromMolarDensity(1., 2., 3., 4., 5.), thickness));
}
}  // namespace

BOOST_AUTO_TEST_CASE(StableKeysSurviveRenumberingAndSubsetLoading) {
  auto a = keyedCylinder(1, "barrel/a");
  auto b = keyedCylinder(2, "barrel/b");
  auto materialA = materialWithThickness(1.);
  auto materialB = materialWithThickness(2.);
  std::vector<const Surface*> sources{a.get(), b.get(), a.get()};
  detail::MaterialSurfaceRegistry registry(sources);
  auto maps = registry.materialMaps(
      {{a->geometryId(), materialA}, {b->geometryId(), materialB}});
  BOOST_CHECK(maps.surfaceMaterials.empty());
  BOOST_CHECK_EQUAL(maps.keyedSurfaces.size(), 2u);
  MaterialMapJsonConverter converter({}, Logging::INFO);
  const auto json = converter.materialMapsToJson(maps);
  auto invalidShape = json;
  invalidShape["KeyedSurfaces"] = nlohmann::json::object();
  BOOST_CHECK_THROW(converter.jsonToMaterialMaps(invalidShape),
                    std::invalid_argument);
  BOOST_CHECK(json.at("KeyedSurfaces").is_array());
  BOOST_CHECK(json.at("Surfaces").at("entries").empty());
  for (const auto& entry : json.at("KeyedSurfaces")) {
    BOOST_CHECK_EQUAL(entry.size(), 3u);
    BOOST_CHECK(entry.contains("key"));
    BOOST_CHECK(entry.contains("geometry_id"));
    BOOST_CHECK(entry.contains("material"));
  }
  TrackingGeometryMaterial loader(converter.jsonToMaterialMaps(json));

  // IDs are swapped, so an ID lookup would silently assign the wrong material.
  auto newA = keyedCylinder(2, "barrel/a");
  auto newB = keyedCylinder(1, "barrel/b");
  std::vector<Surface*> targets{newB.get(), newA.get(), newB.get()};
  loader.apply(targets);
  BOOST_CHECK_EQUAL(
      newA->surfaceMaterial()->materialSlab(Vector2{0., 0.}).thickness(), 1.);
  BOOST_CHECK_EQUAL(
      newB->surfaceMaterial()->materialSlab(Vector2{0., 0.}).thickness(), 2.);
  // Loading consumes the keyed placeholder; real material has no identity API.
  BOOST_CHECK(!detail::materialKey(newA->surfaceMaterial()));

  auto subset = keyedCylinder(17, "barrel/b");
  std::vector<Surface*> subsetSurfaces{subset.get()};
  BOOST_CHECK_NO_THROW(loader.apply(subsetSurfaces));
  BOOST_CHECK_EQUAL(
      subset->surfaceMaterial()->materialSlab(Vector2{0., 0.}).thickness(), 2.);
}

BOOST_AUTO_TEST_CASE(KeyedMapValidationAndLegacyFallback) {
  auto a = keyedCylinder(1, "a");
  auto b = keyedCylinder(2, "b");
  auto material = materialWithThickness(3.);
  TrackingGeometryMaterial maps;
  maps.keyedSurfaces.emplace("a",
                             KeyedSurfaceMaterial{a->geometryId(), material});
  maps.surfaceMaterials.emplace(b->geometryId(), material);
  TrackingGeometryMaterial loader(maps);
  std::vector<Surface*> targets{a.get(), b.get()};
  auto before = a->surfaceMaterialSharedPtr();
  // A missing keyed entry must not fall back to b's matching numeric ID.
  BOOST_CHECK_THROW(loader.apply(targets), std::invalid_argument);
  BOOST_CHECK(a->surfaceMaterialSharedPtr() == before);
  auto duplicate = keyedCylinder(3, "a");
  targets = {a.get(), duplicate.get()};
  BOOST_CHECK_THROW(loader.apply(targets), std::invalid_argument);
  auto duplicateId = keyedCylinder(1, "c");
  targets = {a.get(), duplicateId.get()};
  BOOST_CHECK_THROW(loader.apply(targets), std::invalid_argument);
  // Surface dimensions are not part of the keyed assignment.
  auto resized = keyedCylinder(9, "a", 30.);
  targets = {resized.get()};
  BOOST_CHECK_NO_THROW(loader.apply(targets));
  BOOST_CHECK(resized->surfaceMaterialSharedPtr() == material);

  // Unkeyed targets retain ID-based assignment.
  b->assignSurfaceMaterial(
      std::make_shared<ProtoSurfaceMaterial>(BinUtility{}));
  BOOST_CHECK_NO_THROW(loader.apply(*b));
  BOOST_CHECK(b->surfaceMaterialSharedPtr() == material);

  MaterialMapJsonConverter converter({}, Logging::INFO);
  auto json = converter.materialMapsToJson(maps);
  json["KeyedSurfaces"].push_back(json["KeyedSurfaces"][0]);
  BOOST_CHECK_THROW(converter.jsonToMaterialMaps(json), std::invalid_argument);
  json["KeyedSurfaces"][1]["key"] = "another-geometry/a";
  // Historical numeric IDs can overlap across independent maps.
  BOOST_CHECK_EQUAL(converter.jsonToMaterialMaps(json).keyedSurfaces.size(),
                    2u);
  json["KeyedSurfaces"][1]["key"] = "";
  BOOST_CHECK_THROW(converter.jsonToMaterialMaps(json), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(KeysSurviveProtoAndSurfaceSerialization) {
  auto source = keyedCylinder(1, "stable/key");
  auto payload =
      SurfaceMaterialJsonConverter::toJson(*source->surfaceMaterial());
  auto restored = SurfaceMaterialJsonConverter::fromJson(payload);
  const auto* proto =
      dynamic_cast<const ProtoGridSurfaceMaterial*>(restored.get());
  BOOST_REQUIRE(proto != nullptr);
  BOOST_CHECK_EQUAL(*proto->materialKey(), "stable/key");
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  auto json = SurfaceJsonConverter::toJson(gctx, *source);
  auto surface = SurfaceJsonConverter::fromJson(json);
  BOOST_REQUIRE(detail::materialKey(surface->surfaceMaterial()));
  BOOST_CHECK_EQUAL(*detail::materialKey(surface->surfaceMaterial()),
                    "stable/key");
}

BOOST_AUTO_TEST_CASE(MappingRejectsAmbiguousAndMergedAssignments) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  auto a = keyedCylinder(1, "a");
  auto b = keyedCylinder(2, "a");
  BinnedSurfaceMaterialAccumulator::Config config;
  config.materialSurfaces = {a.get(), b.get()};
  BOOST_CHECK_THROW(BinnedSurfaceMaterialAccumulator(config).createState(gctx),
                    std::invalid_argument);
  b->assignSurfaceMaterial(keyedCylinder(2, "b")->surfaceMaterialSharedPtr());
  b->assignGeometryId(a->geometryId());
  BOOST_CHECK_THROW(BinnedSurfaceMaterialAccumulator(config).createState(gctx),
                    std::invalid_argument);
  b->assignGeometryId(GeometryIdentifier().withSensitive(2));
  config.materialSurfaces.push_back(a.get());
  auto accumulator = std::make_shared<BinnedSurfaceMaterialAccumulator>(config);
  MaterialMapper::Config mapperConfig;
  mapperConfig.surfaceMaterialAccumulator = accumulator;
  mapperConfig.assignmentFinder =
      std::make_shared<IntersectionMaterialAssigner>(
          IntersectionMaterialAssigner::Config{});
  MaterialMapper mapper(mapperConfig);
  auto state = mapper.createState(gctx);
  auto mapped = mapper.finalizeMaps(*state, gctx);
  BOOST_CHECK_EQUAL(mapped.keyedSurfaces.size(), 2u);
  BOOST_CHECK(mapped.surfaceMaterials.empty());
  // Even a no-hit result is an explicit vacuum assignment, not a missing key.
  MaterialMapJsonConverter converter({}, Logging::INFO);
  auto restoredMaps =
      converter.jsonToMaterialMaps(converter.materialMapsToJson(mapped));
  TrackingGeometryMaterial loader(std::move(restoredMaps));
  auto target = keyedCylinder(10, "a");
  BOOST_CHECK_NO_THROW(loader.apply(*target));
  BOOST_CHECK(target->surfaceMaterial()->materialSlab(Vector2{0., 0.}) ==
              MaterialSlab::Nothing());
  b->assignSurfaceMaterial(
      keyedCylinder(2, "changed")->surfaceMaterialSharedPtr());
  BOOST_CHECK_THROW(mapper.finalizeMaps(*state, gctx), std::invalid_argument);
  b->assignSurfaceMaterial(keyedCylinder(2, "b")->surfaceMaterialSharedPtr());

  auto marker = std::make_unique<MergedMaterialMarker>(
      std::vector<MergedMaterialMarker::Origin>{
          {a->geometryId(), "discarded/a"}, {b->geometryId(), "discarded/b"}});
  auto payload = SurfaceMaterialJsonConverter::toJson(*marker);
  auto decoded = SurfaceMaterialJsonConverter::fromJson(payload);
  const auto* decodedMarker =
      dynamic_cast<const MergedMaterialMarker*>(decoded.get());
  BOOST_REQUIRE(decodedMarker);
  BOOST_REQUIRE_EQUAL(decodedMarker->origins().size(), 2u);
  BOOST_CHECK_EQUAL(*decodedMarker->origins()[0].materialKey, "discarded/a");
  BOOST_CHECK(decodedMarker->origins()[1].geometryId == b->geometryId());
  a->assignSurfaceMaterial(std::move(marker));
  BOOST_CHECK_EXCEPTION(
      accumulator->createState(gctx), std::invalid_argument, [](const auto& e) {
        return std::string(e.what()).find("discarded/a") != std::string::npos;
      });
}

BOOST_AUTO_TEST_CASE(KeyedJsonDecoratorRejectsGen1Construction) {
  TrackingGeometryMaterial maps;
  maps.keyedSurfaces.emplace(
      "barrel",
      KeyedSurfaceMaterial{
          {}, std::make_shared<HomogeneousSurfaceMaterial>(MaterialSlab{})});
  MaterialMapJsonConverter converter(MaterialMapJsonConverter::Config{},
                                     Logging::WARNING);
  const std::string path = "keyed-gen1-rejection.json";
  {
    std::ofstream out(path);
    out << converter.materialMapsToJson(maps);
  }
  JsonMaterialDecorator decorator(MaterialMapJsonConverter::Config{}, path,
                                  Logging::WARNING);
  auto world = std::make_shared<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(0., 100., 200.), "world");
  BOOST_CHECK_EXCEPTION(
      TrackingGeometry(world, &decorator), std::invalid_argument,
      [](const auto& error) {
        return std::string(error.what()).find("keyed material map to Gen1") !=
               std::string::npos;
      });
  std::remove(path.c_str());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

ACTS_POP_IGNORE_DEPRECATED()

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/MergedMaterialMarker.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"
#include "ActsPlugins/Json/SurfaceMaterialJsonConverter.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <array>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <variant>
#include <vector>

#include <nlohmann/json.hpp>

using namespace Acts;

namespace ActsTests {

namespace {

/// A few material slabs shared by the index based tests
std::vector<MaterialSlab> testSlabs() {
  std::vector<MaterialSlab> material;
  material.emplace_back(Material::Vacuum(), 0.0);
  material.emplace_back(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0),
                        1.0);
  material.emplace_back(
      Material::fromMolarDensity(11.0, 12.0, 13.0, 14.0, 15.0), 2.0);
  material.emplace_back(
      Material::fromMolarDensity(21.0, 22.0, 23.0, 24.0, 25.0), 3.0);
  return material;
}

/// Grid material is always 2D now, these are the axes the tests bin on
std::unique_ptr<IAxis> testAxis0() {
  return IAxis::createEquidistant(AxisBoundaryType::Bound, -1., 1., 2);
}

std::unique_ptr<IAxis> testAxis1() {
  return IAxis::createEquidistant(AxisBoundaryType::Closed, -std::numbers::pi,
                                  std::numbers::pi, 4);
}

/// Index payload in column major order, i.e. [i0][i1]. Index 0 is vacuum and
/// several bins share an entry, which is what the index storage is for.
std::vector<std::vector<std::size_t>> testIndexPayload() {
  return {std::vector<std::size_t>{1u, 1u, 0u, 2u},
          std::vector<std::size_t>{0u, 3u, 3u, 0u}};
}

std::unique_ptr<GridSurfaceMaterial> makeIndexed() {
  auto axis0 = testAxis0();
  auto axis1 = testAxis1();
  return GridSurfaceMaterial::createIndexed(*axis0, *axis1, testSlabs(),
                                            testIndexPayload());
}

std::unique_ptr<GridSurfaceMaterial> makeGloballyIndexed(
    std::shared_ptr<std::vector<MaterialSlab>> store = nullptr) {
  if (store == nullptr) {
    store = std::make_shared<std::vector<MaterialSlab>>(testSlabs());
  }
  auto axis0 = testAxis0();
  auto axis1 = testAxis1();
  return GridSurfaceMaterial::createGloballyIndexed(
      *axis0, *axis1, std::move(store), testIndexPayload());
}

std::unique_ptr<GridSurfaceMaterial> makeDirect() {
  auto slabs = testSlabs();
  std::vector<std::vector<MaterialSlab>> payload{
      {slabs[1], slabs[1], slabs[0], slabs[2]},
      {slabs[0], slabs[3], slabs[3], slabs[0]}};
  auto axis0 = testAxis0();
  auto axis1 = testAxis1();
  return GridSurfaceMaterial::createDirect(*axis0, *axis1, payload);
}

/// The local points that address the four phi bins of the two z bins
std::vector<Vector2> testPoints() {
  return {{-0.5, -std::numbers::pi * 0.75}, {-0.5, -std::numbers::pi / 4.},
          {-0.5, std::numbers::pi / 4.},    {-0.5, std::numbers::pi * 0.75},
          {0.5, -std::numbers::pi * 0.75},  {0.5, -std::numbers::pi / 4.},
          {0.5, std::numbers::pi / 4.},     {0.5, std::numbers::pi * 0.75}};
}

BinUtility testBinUtility2D() {
  BinUtility bUtility(2, -1., 1., open, AxisDirection::AxisX);
  bUtility += BinUtility(3, -3., 3., open, AxisDirection::AxisY);
  return bUtility;
}

/// The matrix is indexed [bin of the second binning][bin of the first]
MaterialSlabMatrix testMatrix2D() {
  auto slabs = testSlabs();
  MaterialSlabMatrix matrix;
  for (std::size_t i1 = 0; i1 < 3; ++i1) {
    MaterialSlabVector row;
    for (std::size_t i0 = 0; i0 < 2; ++i0) {
      row.push_back(slabs[(i1 * 2 + i0) % slabs.size()]);
    }
    matrix.push_back(std::move(row));
  }
  return matrix;
}

/// Round trip a material through the converter
std::unique_ptr<const ISurfaceMaterial> roundTrip(
    const ISurfaceMaterial& material) {
  nlohmann::json jMaterial = SurfaceMaterialJsonConverter::toJson(material);
  return SurfaceMaterialJsonConverter::fromJson(jMaterial);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(HomogeneousSurfaceMaterialRoundTrip) {
  HomogeneousSurfaceMaterial hsm(
      MaterialSlab(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0), 1.5),
      1., MappingType::PostMapping);

  nlohmann::json jMaterial = SurfaceMaterialJsonConverter::toJson(hsm);
  BOOST_CHECK_EQUAL(jMaterial["type"], "homogeneous");

  auto read = roundTrip(hsm);
  BOOST_REQUIRE(read != nullptr);
  const auto* typed =
      dynamic_cast<const HomogeneousSurfaceMaterial*>(read.get());
  BOOST_REQUIRE(typed != nullptr);
  BOOST_CHECK(typed->mappingType() == MappingType::PostMapping);
  CHECK_CLOSE_ABS(typed->materialSlab(Vector2{0., 0.}).thickness(), 1.5, 1e-5);
}

BOOST_AUTO_TEST_CASE(BinnedSurfaceMaterialRoundTrip) {
  BinnedSurfaceMaterial bsm(testBinUtility2D(), testMatrix2D(), 1.,
                            MappingType::Sensor);

  nlohmann::json jMaterial = SurfaceMaterialJsonConverter::toJson(bsm);
  BOOST_CHECK_EQUAL(jMaterial["type"], "binned");

  auto read = roundTrip(bsm);
  BOOST_REQUIRE(read != nullptr);
  const auto* typed = dynamic_cast<const BinnedSurfaceMaterial*>(read.get());
  BOOST_REQUIRE(typed != nullptr);
  BOOST_CHECK(typed->mappingType() == MappingType::Sensor);
  BOOST_CHECK(typed->binUtility() == bsm.binUtility());

  for (double x : {-0.5, 0.5}) {
    for (double y : {-2., 0., 2.}) {
      Vector2 lp{x, y};
      CHECK_CLOSE_ABS(typed->materialSlab(lp).thickness(),
                      bsm.materialSlab(lp).thickness(), 1e-5);
    }
  }
}

BOOST_AUTO_TEST_CASE(ProtoSurfaceMaterialRoundTrip) {
  ProtoSurfaceMaterial psm(testBinUtility2D(), MappingType::PreMapping);

  nlohmann::json jMaterial = SurfaceMaterialJsonConverter::toJson(psm);
  BOOST_CHECK_EQUAL(jMaterial["type"], "proto");
  BOOST_CHECK_EQUAL(jMaterial["mapMaterial"], true);

  auto read = roundTrip(psm);
  BOOST_REQUIRE(read != nullptr);
  const auto* typed = dynamic_cast<const ProtoSurfaceMaterial*>(read.get());
  BOOST_REQUIRE(typed != nullptr);
  BOOST_CHECK(typed->mappingType() == MappingType::PreMapping);
  BOOST_CHECK(typed->binning() == psm.binning());
}

BOOST_AUTO_TEST_CASE(ProtoGridSurfaceMaterialRoundTrip) {
  MultiAxisSpec2D spec{std::array<AxisSpec, 2u>{
      AxisSpec::Equidistant(4u, -1., 1., AxisBoundaryType::Bound,
                            AxisDirection::AxisX),
      AxisSpec::DeferredVariable({0., 0.25, 1.}, AxisBoundaryType::Bound,
                                 AxisDirection::AxisY)}};
  ProtoGridSurfaceMaterial pgsm(spec, MappingType::PostMapping);

  nlohmann::json jMaterial = SurfaceMaterialJsonConverter::toJson(pgsm);
  BOOST_CHECK_EQUAL(jMaterial["type"], "proto-grid");

  auto read = roundTrip(pgsm);
  BOOST_REQUIRE(read != nullptr);
  const auto* typed = dynamic_cast<const ProtoGridSurfaceMaterial*>(read.get());
  BOOST_REQUIRE(typed != nullptr);
  BOOST_CHECK(typed->mappingType() == MappingType::PostMapping);
  BOOST_CHECK(typed->binning() == pgsm.binning());
}

BOOST_AUTO_TEST_CASE(MergedMaterialMarkerRoundTrip) {
  MergedMaterialMarker marker;

  nlohmann::json jMaterial = SurfaceMaterialJsonConverter::toJson(marker);
  BOOST_CHECK_EQUAL(jMaterial["type"], "merged-material-marker");

  auto read = roundTrip(marker);
  BOOST_REQUIRE(read != nullptr);
  BOOST_CHECK(dynamic_cast<const MergedMaterialMarker*>(read.get()) != nullptr);
}

BOOST_AUTO_TEST_CASE(IndexedGridMaterialRoundTrip) {
  auto ism = makeIndexed();

  nlohmann::json jMaterial = SurfaceMaterialJsonConverter::toJson(*ism);
  BOOST_CHECK_EQUAL(jMaterial["type"], "grid");
  BOOST_CHECK_EQUAL(jMaterial["accessor"]["type"], "indexed");
  BOOST_CHECK_EQUAL(jMaterial["accessor"]["grid"]["axes"].size(), 2u);
  // The delegates are gone with the multi-axis migration
  BOOST_CHECK(!jMaterial.contains("bound_to_grid_local"));
  BOOST_CHECK(!jMaterial.contains("global_to_grid_local"));

  auto read = roundTrip(*ism);
  BOOST_REQUIRE(read != nullptr);
  const auto* typed = dynamic_cast<const GridSurfaceMaterial*>(read.get());
  BOOST_REQUIRE(typed != nullptr);
  const auto& indexed =
      std::get<GridSurfaceMaterial::Indexed>(typed->storage());
  BOOST_CHECK_EQUAL(indexed.material.size(), testSlabs().size());

  for (const Vector2& lp : testPoints()) {
    CHECK_CLOSE_ABS(typed->materialSlab(lp).thickness(),
                    ism->materialSlab(lp).thickness(), 1e-5);
  }
}

BOOST_AUTO_TEST_CASE(GloballyIndexedGridMaterialRoundTrip) {
  auto gism = makeGloballyIndexed();

  nlohmann::json jMaterial = SurfaceMaterialJsonConverter::toJson(*gism);
  BOOST_CHECK_EQUAL(jMaterial["type"], "grid");
  BOOST_CHECK_EQUAL(jMaterial["accessor"]["type"], "globally_indexed");
  // The store travels with the payload, so this reads back standalone. The
  // free-function reader on main refuses this case outright.
  BOOST_CHECK(jMaterial["accessor"].contains("storage_vector"));
  BOOST_CHECK(!jMaterial["accessor"].contains("store"));

  auto read = roundTrip(*gism);
  BOOST_REQUIRE(read != nullptr);
  const auto* typed = dynamic_cast<const GridSurfaceMaterial*>(read.get());
  BOOST_REQUIRE(typed != nullptr);
  // The globally indexed storage must survive, and its store must be filled
  const auto& global =
      std::get<GridSurfaceMaterial::GloballyIndexed>(typed->storage());
  BOOST_REQUIRE(global.material != nullptr);
  BOOST_CHECK_EQUAL(global.material->size(), testSlabs().size());

  for (const Vector2& lp : testPoints()) {
    CHECK_CLOSE_ABS(typed->materialSlab(lp).thickness(),
                    gism->materialSlab(lp).thickness(), 1e-5);
  }
}

BOOST_AUTO_TEST_CASE(GloballyIndexedSharedStoreThroughContext) {
  auto store = std::make_shared<std::vector<MaterialSlab>>(testSlabs());
  auto first = makeGloballyIndexed(store);
  auto second = makeGloballyIndexed(store);

  auto encodeContext =
      SurfaceMaterialJsonConverter::EncodeContext::withStoreTable();
  nlohmann::json jFirst =
      SurfaceMaterialJsonConverter::toJson(*first, encodeContext);
  nlohmann::json jSecond =
      SurfaceMaterialJsonConverter::toJson(*second, encodeContext);

  // One table entry, referenced by both surfaces
  BOOST_REQUIRE_EQUAL(encodeContext.stores().size(), 1u);
  BOOST_CHECK(encodeContext.stores()[0] == store);
  BOOST_CHECK_EQUAL(jFirst["accessor"]["store"], 0u);
  BOOST_CHECK_EQUAL(jSecond["accessor"]["store"], 0u);
  BOOST_CHECK(!jFirst["accessor"].contains("storage_vector"));
  BOOST_CHECK(!jSecond["accessor"].contains("storage_vector"));

  SurfaceMaterialJsonConverter::DecodeContext decodeContext;
  decodeContext.setStores(
      {std::make_shared<std::vector<MaterialSlab>>(testSlabs())});

  auto readFirst =
      SurfaceMaterialJsonConverter::fromJson(jFirst, decodeContext);
  auto readSecond =
      SurfaceMaterialJsonConverter::fromJson(jSecond, decodeContext);
  BOOST_REQUIRE(readFirst != nullptr);
  BOOST_REQUIRE(readSecond != nullptr);

  const auto& globalFirst = std::get<GridSurfaceMaterial::GloballyIndexed>(
      dynamic_cast<const GridSurfaceMaterial&>(*readFirst).storage());
  const auto& globalSecond = std::get<GridSurfaceMaterial::GloballyIndexed>(
      dynamic_cast<const GridSurfaceMaterial&>(*readSecond).storage());
  // The sharing must survive the round trip
  BOOST_CHECK(globalFirst.material == globalSecond.material);
  BOOST_CHECK(globalFirst.material == decodeContext.store(0u));
}

BOOST_AUTO_TEST_CASE(DistinctStoresGetSequentialIds) {
  auto storeA = std::make_shared<std::vector<MaterialSlab>>(testSlabs());
  // Same content, different allocation: stores are keyed on identity
  auto storeB = std::make_shared<std::vector<MaterialSlab>>(testSlabs());

  auto ctx = SurfaceMaterialJsonConverter::EncodeContext::withStoreTable();
  nlohmann::json jA =
      SurfaceMaterialJsonConverter::toJson(*makeGloballyIndexed(storeA), ctx);
  nlohmann::json jB =
      SurfaceMaterialJsonConverter::toJson(*makeGloballyIndexed(storeB), ctx);
  nlohmann::json jA2 =
      SurfaceMaterialJsonConverter::toJson(*makeGloballyIndexed(storeA), ctx);

  BOOST_CHECK_EQUAL(jA["accessor"]["store"], 0u);
  BOOST_CHECK_EQUAL(jB["accessor"]["store"], 1u);
  BOOST_CHECK_EQUAL(jA2["accessor"]["store"], 0u);
  BOOST_REQUIRE_EQUAL(ctx.stores().size(), 2u);
  BOOST_CHECK(ctx.stores()[0] == storeA);
  BOOST_CHECK(ctx.stores()[1] == storeB);
}

BOOST_AUTO_TEST_CASE(StoreReferenceWithoutTableThrows) {
  auto store = std::make_shared<std::vector<MaterialSlab>>(testSlabs());
  auto ctx = SurfaceMaterialJsonConverter::EncodeContext::withStoreTable();
  nlohmann::json jMaterial =
      SurfaceMaterialJsonConverter::toJson(*makeGloballyIndexed(store), ctx);

  // A default decode context has no table at all
  BOOST_CHECK_THROW(SurfaceMaterialJsonConverter::fromJson(jMaterial),
                    std::invalid_argument);

  // A table that does not reach the referenced id
  SurfaceMaterialJsonConverter::DecodeContext empty;
  empty.setStores({});
  BOOST_CHECK_THROW(SurfaceMaterialJsonConverter::fromJson(jMaterial, empty),
                    std::invalid_argument);

  nlohmann::json jOutOfRange = jMaterial;
  jOutOfRange["accessor"]["store"] = 7u;
  SurfaceMaterialJsonConverter::DecodeContext oneEntry;
  oneEntry.setStores({store});
  BOOST_CHECK_THROW(
      SurfaceMaterialJsonConverter::fromJson(jOutOfRange, oneEntry),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(DirectGridMaterialRoundTrip) {
  auto gsm = makeDirect();

  nlohmann::json jMaterial = SurfaceMaterialJsonConverter::toJson(*gsm);
  BOOST_CHECK_EQUAL(jMaterial["type"], "grid");
  BOOST_CHECK_EQUAL(jMaterial["accessor"]["type"], "direct");
  BOOST_CHECK(!jMaterial["accessor"].contains("storage_vector"));

  auto read = roundTrip(*gsm);
  BOOST_REQUIRE(read != nullptr);
  const auto* typed = dynamic_cast<const GridSurfaceMaterial*>(read.get());
  BOOST_REQUIRE(typed != nullptr);
  BOOST_CHECK(
      std::holds_alternative<GridSurfaceMaterial::Direct>(typed->storage()));

  for (const Vector2& lp : testPoints()) {
    CHECK_CLOSE_ABS(typed->materialSlab(lp).thickness(),
                    gsm->materialSlab(lp).thickness(), 1e-5);
  }
}

BOOST_AUTO_TEST_CASE(EncoderCoversAllSurfaceMaterials) {
  const auto& cfg = SurfaceMaterialJsonConverter::defaultConfig();

  std::vector<std::shared_ptr<const ISurfaceMaterial>> materials;
  materials.push_back(std::make_shared<const HomogeneousSurfaceMaterial>(
      MaterialSlab(Material::fromMolarDensity(1.0, 2.0, 3.0, 4.0, 5.0), 1.)));
  materials.push_back(std::make_shared<const BinnedSurfaceMaterial>(
      testBinUtility2D(), testMatrix2D()));
  materials.push_back(
      std::make_shared<const ProtoSurfaceMaterial>(testBinUtility2D()));
  materials.push_back(std::make_shared<const ProtoGridSurfaceMaterial>(
      MultiAxisSpec2D{std::array<AxisSpec, 2u>{
          AxisSpec::Equidistant(2u, 0., 1., AxisBoundaryType::Bound,
                                AxisDirection::AxisX),
          AxisSpec::Equidistant(2u, 0., 1., AxisBoundaryType::Bound,
                                AxisDirection::AxisY)}}));
  materials.push_back(std::make_shared<const MergedMaterialMarker>());
  materials.push_back(makeIndexed());
  materials.push_back(makeGloballyIndexed());
  materials.push_back(makeDirect());

  for (std::size_t im = 0; im < materials.size(); ++im) {
    const auto& material = materials[im];
    BOOST_TEST_CONTEXT("material " << im) {
      // Exactly one encoder must claim the type, otherwise the call throws
      BOOST_CHECK(cfg.encoder.hasFunction(*material));
      nlohmann::json jMaterial;
      BOOST_REQUIRE_NO_THROW(
          jMaterial = SurfaceMaterialJsonConverter::toJson(*material));
      // Every tag the encoder can emit must be known to the decoder
      BOOST_CHECK(cfg.decoder.hasKind(jMaterial["type"].get<std::string>()));
    }
  }

  // One concrete class carries the whole grid material family
  BOOST_CHECK(cfg.encoder.hasFunction<HomogeneousSurfaceMaterial>());
  BOOST_CHECK(cfg.encoder.hasFunction<BinnedSurfaceMaterial>());
  BOOST_CHECK(cfg.encoder.hasFunction<ProtoSurfaceMaterial>());
  BOOST_CHECK(cfg.encoder.hasFunction<ProtoGridSurfaceMaterial>());
  BOOST_CHECK(cfg.encoder.hasFunction<MergedMaterialMarker>());
  BOOST_CHECK(cfg.encoder.hasFunction<GridSurfaceMaterial>());
  BOOST_CHECK_EQUAL(cfg.encoder.size(), 6u);
  BOOST_CHECK_EQUAL(cfg.decoder.size(), 6u);
}

BOOST_AUTO_TEST_CASE(MissingAndUnknownTypeTagThrow) {
  nlohmann::json jMissing;
  jMissing["mapMaterial"] = true;
  BOOST_CHECK_THROW(SurfaceMaterialJsonConverter::fromJson(jMissing),
                    std::invalid_argument);

  nlohmann::json jUnknown;
  jUnknown["mapMaterial"] = true;
  jUnknown["type"] = "not-a-material";
  BOOST_CHECK_THROW(SurfaceMaterialJsonConverter::fromJson(jUnknown),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(UnmappedMaterialYieldsNoMaterial) {
  nlohmann::json jMaterial;
  jMaterial["type"] = "proto";
  jMaterial["mapMaterial"] = false;
  BOOST_CHECK(SurfaceMaterialJsonConverter::fromJson(jMaterial) == nullptr);

  // A proto material without binning is flagged as not mapped on write
  ProtoSurfaceMaterial psm{BinUtility{}};
  nlohmann::json jProto = SurfaceMaterialJsonConverter::toJson(psm);
  BOOST_CHECK_EQUAL(jProto["mapMaterial"], false);
  BOOST_CHECK(SurfaceMaterialJsonConverter::fromJson(jProto) == nullptr);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

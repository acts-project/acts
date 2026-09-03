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
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"
#include "ActsPlugins/Json/SurfaceMaterialJsonConverter.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <array>
#include <memory>
#include <numbers>
#include <stdexcept>
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

using EqBound = GridAxisGenerators::EqBound;
using EqBoundEqClosed = GridAxisGenerators::EqBoundEqClosed;

using IndexEqGrid = EqBound::grid_type<std::size_t>;
using IndexEqEqGrid = EqBoundEqClosed::grid_type<std::size_t>;
using SlabEqGrid = EqBound::grid_type<MaterialSlab>;

template <typename material_t>
typename material_t::BoundToGridLocalDelegate boundToLocal1D() {
  typename material_t::BoundToGridLocalDelegate delegate;
  delegate.template connect<&GridAccess::LocalSubspace<0u>::toGridLocal>(
      std::make_unique<const GridAccess::LocalSubspace<0u>>());
  return delegate;
}

template <typename material_t>
typename material_t::GlobalToGridLocalDelegate globalToLocal1D() {
  typename material_t::GlobalToGridLocalDelegate delegate;
  delegate.template connect<
      &GridAccess::GlobalSubspace<AxisDirection::AxisX>::toGridLocal>(
      std::make_unique<
          const GridAccess::GlobalSubspace<AxisDirection::AxisX>>());
  return delegate;
}

template <typename material_t>
typename material_t::BoundToGridLocalDelegate boundToLocal2D() {
  typename material_t::BoundToGridLocalDelegate delegate;
  delegate.template connect<&GridAccess::LocalSubspace<0u, 1u>::toGridLocal>(
      std::make_unique<const GridAccess::LocalSubspace<0u, 1u>>());
  return delegate;
}

template <typename material_t>
typename material_t::GlobalToGridLocalDelegate globalToLocal2D() {
  using Subspace =
      GridAccess::GlobalSubspace<AxisDirection::AxisZ, AxisDirection::AxisPhi>;
  typename material_t::GlobalToGridLocalDelegate delegate;
  delegate.template connect<&Subspace::toGridLocal>(
      std::make_unique<const Subspace>());
  return delegate;
}

/// Fill a 1D index grid with the entries 1, 0, 2, 2, 3
IndexEqGrid filledIndexGrid1D() {
  EqBound eqBound{{0., 5.}, 5};
  IndexEqGrid grid{eqBound()};
  grid.atPosition(IndexEqGrid::point_t{0.5}) = 1u;
  grid.atPosition(IndexEqGrid::point_t{1.5}) = 0u;
  grid.atPosition(IndexEqGrid::point_t{2.5}) = 2u;
  grid.atPosition(IndexEqGrid::point_t{3.5}) = 2u;
  grid.atPosition(IndexEqGrid::point_t{4.5}) = 3u;
  return grid;
}

std::shared_ptr<const IndexedSurfaceMaterial<IndexEqGrid>> makeIndexed1D() {
  using Material_t = IndexedSurfaceMaterial<IndexEqGrid>;
  return std::make_shared<const Material_t>(
      filledIndexGrid1D(), IndexedMaterialAccessor{testSlabs()},
      boundToLocal1D<Material_t>(), globalToLocal1D<Material_t>());
}

std::shared_ptr<const GloballyIndexedSurfaceMaterial<IndexEqGrid>>
makeGloballyIndexed1D(MaterialSlabStore store = nullptr, bool shared = false) {
  using Material_t = GloballyIndexedSurfaceMaterial<IndexEqGrid>;
  if (store == nullptr) {
    store = std::make_shared<std::vector<MaterialSlab>>(testSlabs());
  }
  return std::make_shared<const Material_t>(
      filledIndexGrid1D(),
      GloballyIndexedMaterialAccessor{std::move(store), shared},
      boundToLocal1D<Material_t>(), globalToLocal1D<Material_t>());
}

std::shared_ptr<const IndexedSurfaceMaterial<IndexEqEqGrid>> makeIndexed2D() {
  using Material_t = IndexedSurfaceMaterial<IndexEqEqGrid>;
  using Point = IndexEqEqGrid::point_t;

  EqBoundEqClosed eqeqBound{
      {-1., 1.}, 2, {-std::numbers::pi, std::numbers::pi}, 4};
  IndexEqEqGrid grid{eqeqBound()};
  grid.atPosition(Point{-0.5, -std::numbers::pi * 0.75}) = 1u;
  grid.atPosition(Point{-0.5, -std::numbers::pi / 4.}) = 1u;
  grid.atPosition(Point{-0.5, std::numbers::pi / 4.}) = 0u;
  grid.atPosition(Point{-0.5, std::numbers::pi * 0.75}) = 2u;
  grid.atPosition(Point{0.5, -std::numbers::pi * 0.75}) = 0u;
  grid.atPosition(Point{0.5, -std::numbers::pi / 4.}) = 3u;
  grid.atPosition(Point{0.5, std::numbers::pi / 4.}) = 3u;
  grid.atPosition(Point{0.5, std::numbers::pi * 0.75}) = 0u;

  return std::make_shared<const Material_t>(
      std::move(grid), IndexedMaterialAccessor{testSlabs()},
      boundToLocal2D<Material_t>(), globalToLocal2D<Material_t>());
}

std::shared_ptr<const GridSurfaceMaterial<SlabEqGrid>> makeSlabGrid1D() {
  using Material_t = GridSurfaceMaterial<SlabEqGrid>;
  auto slabs = testSlabs();

  EqBound eqBound{{0., 4.}, 4};
  SlabEqGrid grid{eqBound()};
  for (std::size_t ib = 0; ib < 4; ++ib) {
    grid.atPosition(SlabEqGrid::point_t{ib + 0.5}) = slabs[ib];
  }
  return std::make_shared<const Material_t>(
      std::move(grid), GridMaterialAccessor{}, boundToLocal1D<Material_t>(),
      globalToLocal1D<Material_t>());
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

BOOST_AUTO_TEST_CASE(IndexedSurfaceMaterial1DRoundTrip) {
  auto ism = makeIndexed1D();

  nlohmann::json jMaterial = SurfaceMaterialJsonConverter::toJson(*ism);
  BOOST_CHECK_EQUAL(jMaterial["type"], "grid");
  BOOST_CHECK_EQUAL(jMaterial["accessor"]["type"], "indexed");

  auto read = roundTrip(*ism);
  BOOST_REQUIRE(read != nullptr);
  const auto* typed =
      dynamic_cast<const IGridSurfaceMaterial<std::size_t>*>(read.get());
  BOOST_REQUIRE(typed != nullptr);
  BOOST_CHECK(dynamic_cast<const IndexedMaterialAccessor*>(
                  &typed->materialAccessor()) != nullptr);

  for (double x : {0.5, 1.5, 2.5, 3.5, 4.5}) {
    Vector2 lp{x, 0.};
    CHECK_CLOSE_ABS(read->materialSlab(lp).thickness(),
                    ism->materialSlab(lp).thickness(), 1e-5);
  }
}

BOOST_AUTO_TEST_CASE(IndexedSurfaceMaterial2DRoundTrip) {
  auto ism = makeIndexed2D();

  nlohmann::json jMaterial = SurfaceMaterialJsonConverter::toJson(*ism);
  BOOST_CHECK_EQUAL(jMaterial["type"], "grid");
  BOOST_CHECK_EQUAL(jMaterial["accessor"]["type"], "indexed");
  BOOST_CHECK_EQUAL(jMaterial["accessor"]["grid"]["axes"].size(), 2u);

  auto read = roundTrip(*ism);
  BOOST_REQUIRE(read != nullptr);

  for (double z : {-0.5, 0.5}) {
    for (double phi : {-std::numbers::pi * 0.75, -std::numbers::pi / 4.,
                       std::numbers::pi / 4., std::numbers::pi * 0.75}) {
      Vector2 lp{z, phi};
      CHECK_CLOSE_ABS(read->materialSlab(lp).thickness(),
                      ism->materialSlab(lp).thickness(), 1e-5);
    }
  }
}

BOOST_AUTO_TEST_CASE(GloballyIndexedSurfaceMaterialRoundTrip) {
  auto gism = makeGloballyIndexed1D(nullptr, true);

  nlohmann::json jMaterial = SurfaceMaterialJsonConverter::toJson(*gism);
  BOOST_CHECK_EQUAL(jMaterial["type"], "grid");
  BOOST_CHECK_EQUAL(jMaterial["accessor"]["type"], "globally_indexed");
  // Without a store table the payload is self-contained
  BOOST_CHECK(jMaterial["accessor"].contains("storage_vector"));
  BOOST_CHECK(!jMaterial["accessor"].contains("store"));
  BOOST_CHECK_EQUAL(jMaterial["accessor"]["shared_entries"], true);

  auto read = roundTrip(*gism);
  BOOST_REQUIRE(read != nullptr);
  const auto* typed =
      dynamic_cast<const IGridSurfaceMaterial<std::size_t>*>(read.get());
  BOOST_REQUIRE(typed != nullptr);
  // The globally indexed accessor must survive, it used to be read back as a
  // locally indexed one
  const auto* accessor = dynamic_cast<const GloballyIndexedMaterialAccessor*>(
      &typed->materialAccessor());
  BOOST_REQUIRE(accessor != nullptr);
  BOOST_CHECK(accessor->sharedEntries);
  // The inlined store used to come back empty
  BOOST_REQUIRE(accessor->slabStore != nullptr);
  BOOST_CHECK_EQUAL(accessor->slabStore->size(), testSlabs().size());

  auto gridView = typed->gridConstView();
  for (std::size_t ib = 1; ib <= 5; ++ib) {
    BOOST_CHECK_EQUAL(gridView.atLocalBins({ib}),
                      gism->gridConstView().atLocalBins({ib}));
  }

  for (double x : {0.5, 1.5, 2.5, 3.5, 4.5}) {
    Vector2 lp{x, 0.};
    CHECK_CLOSE_ABS(read->materialSlab(lp).thickness(),
                    gism->materialSlab(lp).thickness(), 1e-5);
  }
}

BOOST_AUTO_TEST_CASE(GloballyIndexedSharedStoreThroughContext) {
  auto store = std::make_shared<std::vector<MaterialSlab>>(testSlabs());
  auto first = makeGloballyIndexed1D(store);
  auto second = makeGloballyIndexed1D(store);

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

  const auto* accessorFirst =
      dynamic_cast<const GloballyIndexedMaterialAccessor*>(
          &dynamic_cast<const IGridSurfaceMaterial<std::size_t>&>(*readFirst)
               .materialAccessor());
  const auto* accessorSecond =
      dynamic_cast<const GloballyIndexedMaterialAccessor*>(
          &dynamic_cast<const IGridSurfaceMaterial<std::size_t>&>(*readSecond)
               .materialAccessor());
  BOOST_REQUIRE(accessorFirst != nullptr);
  BOOST_REQUIRE(accessorSecond != nullptr);
  // The sharing must survive the round trip
  BOOST_CHECK(accessorFirst->slabStore == accessorSecond->slabStore);
  BOOST_CHECK(accessorFirst->slabStore == decodeContext.store(0u));
}

BOOST_AUTO_TEST_CASE(DistinctStoresGetSequentialIds) {
  auto storeA = std::make_shared<std::vector<MaterialSlab>>(testSlabs());
  // Same content, different allocation: stores are keyed on identity
  auto storeB = std::make_shared<std::vector<MaterialSlab>>(testSlabs());

  auto ctx = SurfaceMaterialJsonConverter::EncodeContext::withStoreTable();
  nlohmann::json jA =
      SurfaceMaterialJsonConverter::toJson(*makeGloballyIndexed1D(storeA), ctx);
  nlohmann::json jB =
      SurfaceMaterialJsonConverter::toJson(*makeGloballyIndexed1D(storeB), ctx);
  nlohmann::json jA2 =
      SurfaceMaterialJsonConverter::toJson(*makeGloballyIndexed1D(storeA), ctx);

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
      SurfaceMaterialJsonConverter::toJson(*makeGloballyIndexed1D(store), ctx);

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

BOOST_AUTO_TEST_CASE(GridSurfaceMaterialRoundTrip) {
  auto gsm = makeSlabGrid1D();

  nlohmann::json jMaterial = SurfaceMaterialJsonConverter::toJson(*gsm);
  BOOST_CHECK_EQUAL(jMaterial["type"], "grid");
  BOOST_CHECK_EQUAL(jMaterial["accessor"]["type"], "grid");

  auto read = roundTrip(*gsm);
  BOOST_REQUIRE(read != nullptr);
  BOOST_CHECK(dynamic_cast<const IGridSurfaceMaterial<MaterialSlab>*>(
                  read.get()) != nullptr);

  for (double x : {0.5, 1.5, 2.5, 3.5}) {
    Vector2 lp{x, 0.};
    CHECK_CLOSE_ABS(read->materialSlab(lp).thickness(),
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
  materials.push_back(makeIndexed1D());
  materials.push_back(makeIndexed2D());
  materials.push_back(makeGloballyIndexed1D());
  materials.push_back(makeSlabGrid1D());

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

  // The abstract grid bases carry the whole grid material family
  BOOST_CHECK(cfg.encoder.hasFunction<HomogeneousSurfaceMaterial>());
  BOOST_CHECK(cfg.encoder.hasFunction<BinnedSurfaceMaterial>());
  BOOST_CHECK(cfg.encoder.hasFunction<ProtoSurfaceMaterial>());
  BOOST_CHECK(cfg.encoder.hasFunction<ProtoGridSurfaceMaterial>());
  BOOST_CHECK(cfg.encoder.hasFunction<MergedMaterialMarker>());
  BOOST_CHECK(cfg.encoder.hasFunction<IGridSurfaceMaterial<MaterialSlab>>());
  BOOST_CHECK(cfg.encoder.hasFunction<IGridSurfaceMaterial<std::size_t>>());
  BOOST_CHECK_EQUAL(cfg.encoder.size(), 7u);
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

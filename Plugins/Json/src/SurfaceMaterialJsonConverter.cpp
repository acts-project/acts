// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/SurfaceMaterialJsonConverter.hpp"

#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterialFactory.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/MergedMaterialMarker.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Utilities/AnyGridView.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/IGrid.hpp"
#include "ActsPlugins/Json/AxisSpecJsonConverter.hpp"
#include "ActsPlugins/Json/GridJsonConverter.hpp"
#include "ActsPlugins/Json/MaterialJsonConverter.hpp"
#include "ActsPlugins/Json/UtilitiesJsonConverter.hpp"

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace Acts;

using EncodeContext = SurfaceMaterialJsonConverter::EncodeContext;
using DecodeContext = SurfaceMaterialJsonConverter::DecodeContext;

/// Payload type tags, shared between the encoder and the decoder
constexpr const char* kHomogeneousTag = "homogeneous";
constexpr const char* kBinnedTag = "binned";
constexpr const char* kProtoTag = "proto";
constexpr const char* kProtoGridTag = "proto-grid";
constexpr const char* kMergedMarkerTag = "merged-material-marker";
constexpr const char* kGridTag = "grid";

/// Grid accessor tags
constexpr const char* kGridAccessorTag = "grid";
constexpr const char* kIndexedAccessorTag = "indexed";
constexpr const char* kGloballyIndexedAccessorTag = "globally_indexed";

nlohmann::json homogeneousToJson(const HomogeneousSurfaceMaterial& material,
                                 EncodeContext& /*ctx*/) {
  nlohmann::json jMaterial;
  jMaterial[jsonKey().typekey] = kHomogeneousTag;
  jMaterial[jsonKey().maptype] = nlohmann::json(material.mappingType());
  jMaterial[jsonKey().mapkey] = true;
  nlohmann::json jSlab(material.materialSlab());
  jMaterial[jsonKey().datakey] = nlohmann::json::array({
      nlohmann::json::array({
          jSlab,
      }),
  });
  return jMaterial;
}

nlohmann::json binnedToJson(const BinnedSurfaceMaterial& material,
                            EncodeContext& /*ctx*/) {
  nlohmann::json jMaterial;
  jMaterial[jsonKey().typekey] = kBinnedTag;
  jMaterial[jsonKey().maptype] = nlohmann::json(material.mappingType());
  jMaterial[jsonKey().mapkey] = true;

  nlohmann::json jMatrix = nlohmann::json::array();
  for (const auto& mpVector : material.fullMaterial()) {
    nlohmann::json jVector = nlohmann::json::array();
    for (const auto& mp : mpVector) {
      jVector.push_back(nlohmann::json(mp));
    }
    jMatrix.push_back(std::move(jVector));
  }
  jMaterial[jsonKey().datakey] = std::move(jMatrix);
  jMaterial[jsonKey().binkey] = nlohmann::json(material.binUtility());
  return jMaterial;
}

nlohmann::json protoToJson(const ProtoSurfaceMaterial& material,
                           EncodeContext& /*ctx*/) {
  nlohmann::json jMaterial;
  jMaterial[jsonKey().typekey] = kProtoTag;
  jMaterial[jsonKey().maptype] = nlohmann::json(material.mappingType());
  // A proto material without any actual binning is not mapped onto
  jMaterial[jsonKey().mapkey] = false;
  const BinUtility& bUtility = material.binning();
  for (const auto& bData : bUtility.binningData()) {
    if (bData.bins() > 1) {
      jMaterial[jsonKey().mapkey] = true;
      break;
    }
  }
  jMaterial[jsonKey().binkey] = nlohmann::json(bUtility);
  return jMaterial;
}

nlohmann::json protoGridToJson(const ProtoGridSurfaceMaterial& material,
                               EncodeContext& /*ctx*/) {
  nlohmann::json jMaterial;
  jMaterial[jsonKey().typekey] = kProtoGridTag;
  jMaterial[jsonKey().maptype] = nlohmann::json(material.mappingType());
  jMaterial[jsonKey().mapkey] = true;
  jMaterial["axis_specs"] =
      MultiAxisSpecJsonConverter::toJson(material.binning());
  return jMaterial;
}

nlohmann::json mergedMarkerToJson(const MergedMaterialMarker& /*material*/,
                                  EncodeContext& /*ctx*/) {
  nlohmann::json jMaterial;
  jMaterial[jsonKey().typekey] = kMergedMarkerTag;
  // Flag as "mapped" so the reader does not discard it
  jMaterial[jsonKey().mapkey] = true;
  return jMaterial;
}

/// Type-erased grid serialization, mirroring GridJsonConverter::toJson
template <typename value_t>
nlohmann::json gridToJson(const IGrid& grid,
                          const AnyGridConstView<value_t>& view) {
  nlohmann::json jGrid;

  auto axes = grid.axes();
  nlohmann::json jAxes = nlohmann::json::array();
  for (const IAxis* axis : axes) {
    jAxes.push_back(AxisJsonConverter::toJson(*axis));
  }
  jGrid["axes"] = jAxes;

  nlohmann::json jData = nlohmann::json::array();
  if (axes.size() == 1u) {
    for (std::size_t ib0 = 1u; ib0 <= axes[0u]->getNBins(); ++ib0) {
      nlohmann::json jEntry = nlohmann::json::array();
      jEntry.push_back(std::array<std::size_t, 1u>{ib0});
      jEntry.push_back(nlohmann::json(view.atLocalBins({ib0})));
      jData.push_back(std::move(jEntry));
    }
  } else if (axes.size() == 2u) {
    for (std::size_t ib0 = 1u; ib0 <= axes[0u]->getNBins(); ++ib0) {
      for (std::size_t ib1 = 1u; ib1 <= axes[1u]->getNBins(); ++ib1) {
        nlohmann::json jEntry = nlohmann::json::array();
        jEntry.push_back(std::array<std::size_t, 2u>{ib0, ib1});
        jEntry.push_back(nlohmann::json(view.atLocalBins({ib0, ib1})));
        jData.push_back(std::move(jEntry));
      }
    }
  } else {
    throw std::invalid_argument(
        "SurfaceMaterialJsonConverter: only 1D and 2D material grids are "
        "supported");
  }
  jGrid["data"] = jData;
  return jGrid;
}

nlohmann::json slabsToJson(const std::vector<MaterialSlab>& slabs) {
  nlohmann::json jSlabs = nlohmann::json::array();
  for (const auto& msl : slabs) {
    jSlabs.push_back(nlohmann::json(msl));
  }
  return jSlabs;
}

/// Encoder for the whole grid material family: the concrete grid shape and
/// the accessor type are resolved at runtime through the abstract base.
template <typename value_t>
nlohmann::json gridMaterialToJson(const IGridSurfaceMaterial<value_t>& material,
                                  const std::string& accessorType,
                                  EncodeContext& ctx) {
  nlohmann::json jMaterial;
  jMaterial[jsonKey().typekey] = kGridTag;
  jMaterial[jsonKey().mapkey] = true;

  nlohmann::json jAccessor;
  jAccessor["type"] = accessorType;
  if (accessorType == kIndexedAccessorTag) {
    const auto& accessor = dynamic_cast<const IndexedMaterialAccessor&>(
        material.materialAccessor());
    jAccessor["storage_vector"] = slabsToJson(accessor.material);
  } else if (accessorType == kGloballyIndexedAccessorTag) {
    const auto& accessor = dynamic_cast<const GloballyIndexedMaterialAccessor&>(
        material.materialAccessor());
    if (accessor.slabStore == nullptr) {
      throw std::invalid_argument(
          "SurfaceMaterialJsonConverter: globally indexed material without a "
          "slab store");
    }
    jAccessor["shared_entries"] = accessor.sharedEntries;
    if (ctx.storeTableEnabled()) {
      // The store lives once in the document, the entry only references it
      jAccessor["store"] = ctx.storeId(accessor.slabStore);
    } else {
      // Standalone payload, inline the store to keep it self-contained
      jAccessor["storage_vector"] = slabsToJson(*accessor.slabStore);
    }
  }
  jAccessor["grid"] =
      gridToJson<value_t>(material.grid(), material.gridConstView());
  jMaterial["accessor"] = std::move(jAccessor);

  jMaterial["global_to_grid_local"] =
      GridAccessJsonConverter::toJson(material.globalToGridLocal());
  jMaterial["bound_to_grid_local"] =
      GridAccessJsonConverter::toJson(material.boundToGridLocal());
  return jMaterial;
}

nlohmann::json slabGridMaterialToJson(
    const IGridSurfaceMaterial<MaterialSlab>& material, EncodeContext& ctx) {
  return gridMaterialToJson<MaterialSlab>(material, kGridAccessorTag, ctx);
}

nlohmann::json indexGridMaterialToJson(
    const IGridSurfaceMaterial<std::size_t>& material, EncodeContext& ctx) {
  // Both index-based accessors share the grid value type, they are told
  // apart by the accessor instance
  const IGridMaterialAccessor& accessor = material.materialAccessor();
  if (dynamic_cast<const IndexedMaterialAccessor*>(&accessor) != nullptr) {
    return gridMaterialToJson<std::size_t>(material, kIndexedAccessorTag, ctx);
  }
  if (dynamic_cast<const GloballyIndexedMaterialAccessor*>(&accessor) !=
      nullptr) {
    return gridMaterialToJson<std::size_t>(material,
                                           kGloballyIndexedAccessorTag, ctx);
  }
  throw std::invalid_argument(
      "SurfaceMaterialJsonConverter: unknown index grid material accessor");
}

MappingType readMappingType(const nlohmann::json& jMaterial) {
  MappingType mappingType = MappingType::Default;
  if (jMaterial.contains(jsonKey().maptype) &&
      !jMaterial.at(jsonKey().maptype).is_null()) {
    from_json(jMaterial.at(jsonKey().maptype), mappingType);
  }
  return mappingType;
}

std::unique_ptr<const ISurfaceMaterial> homogeneousFromJson(
    const nlohmann::json& jMaterial, const DecodeContext& /*ctx*/) {
  MaterialSlabMatrix matrix;
  from_json(jMaterial.at(jsonKey().datakey), matrix);
  if (matrix.empty() || matrix[0].empty()) {
    throw std::invalid_argument(
        "SurfaceMaterialJsonConverter: homogeneous material without data");
  }
  return std::make_unique<const HomogeneousSurfaceMaterial>(
      matrix[0][0], 1., readMappingType(jMaterial));
}

std::unique_ptr<const ISurfaceMaterial> binnedFromJson(
    const nlohmann::json& jMaterial, const DecodeContext& /*ctx*/) {
  BinUtility bUtility;
  from_json(jMaterial.at(jsonKey().binkey), bUtility);
  MaterialSlabMatrix matrix;
  from_json(jMaterial.at(jsonKey().datakey), matrix);
  return std::make_unique<const BinnedSurfaceMaterial>(
      bUtility, std::move(matrix), 1., readMappingType(jMaterial));
}

std::unique_ptr<const ISurfaceMaterial> protoFromJson(
    const nlohmann::json& jMaterial, const DecodeContext& /*ctx*/) {
  BinUtility bUtility;
  if (jMaterial.contains(jsonKey().binkey) &&
      !jMaterial.at(jsonKey().binkey).is_null()) {
    from_json(jMaterial.at(jsonKey().binkey), bUtility);
  }
  return std::make_unique<const ProtoSurfaceMaterial>(
      bUtility, readMappingType(jMaterial));
}

std::unique_ptr<const ISurfaceMaterial> protoGridFromJson(
    const nlohmann::json& jMaterial, const DecodeContext& /*ctx*/) {
  MultiAxisSpec spec =
      MultiAxisSpecJsonConverter::fromJson(jMaterial.at("axis_specs"));
  if (spec.size() != 2u) {
    throw std::invalid_argument(
        "SurfaceMaterialJsonConverter: proto grid material needs exactly two "
        "axis specs");
  }
  MultiAxisSpec2D spec2D{
      std::array<AxisSpec, 2u>{spec.axisSpec(0u), spec.axisSpec(1u)}};
  return std::make_unique<const ProtoGridSurfaceMaterial>(
      spec2D, readMappingType(jMaterial));
}

std::unique_ptr<const ISurfaceMaterial> mergedMarkerFromJson(
    const nlohmann::json& /*jMaterial*/, const DecodeContext& /*ctx*/) {
  return std::make_unique<const MergedMaterialMarker>();
}

/// Read the 1D grid payload in axis order, the json holds 1-based local bins
template <typename value_t>
std::vector<value_t> gridPayload1D(const nlohmann::json& jGrid,
                                   std::size_t nBins, const value_t& fill) {
  std::vector<value_t> payload(nBins, fill);
  for (const auto& jEntry : jGrid.at("data")) {
    std::array<std::size_t, 1u> lBin = jEntry.at(0u);
    if (jEntry.at(1u).is_null() || lBin[0u] < 1u || lBin[0u] > nBins) {
      continue;
    }
    payload[lBin[0u] - 1u] = jEntry.at(1u).template get<value_t>();
  }
  return payload;
}

/// Read the 2D grid payload in column major order, i.e. [i0][i1]
template <typename value_t>
std::vector<std::vector<value_t>> gridPayload2D(const nlohmann::json& jGrid,
                                                std::size_t nBins0,
                                                std::size_t nBins1,
                                                const value_t& fill) {
  std::vector<std::vector<value_t>> payload(nBins0,
                                            std::vector<value_t>(nBins1, fill));
  for (const auto& jEntry : jGrid.at("data")) {
    std::array<std::size_t, 2u> lBin = jEntry.at(0u);
    if (jEntry.at(1u).is_null() || lBin[0u] < 1u || lBin[0u] > nBins0 ||
        lBin[1u] < 1u || lBin[1u] > nBins1) {
      continue;
    }
    payload[lBin[0u] - 1u][lBin[1u] - 1u] =
        jEntry.at(1u).template get<value_t>();
  }
  return payload;
}

/// Build the grid material from the type-erased axes and the given accessor
template <typename accessor_t>
std::unique_ptr<const ISurfaceMaterial> createGridMaterial(
    const nlohmann::json& jMaterial, const nlohmann::json& jGrid,
    accessor_t&& accessor, const typename accessor_t::grid_value_type& fill) {
  const nlohmann::json& jAxes = jGrid.at("axes");
  if (jAxes.size() == 1u) {
    std::unique_ptr<IAxis> axis = AxisJsonConverter::fromJson(jAxes.at(0u));
    auto payload = gridPayload1D<typename accessor_t::grid_value_type>(
        jGrid, axis->getNBins(), fill);
    return GridSurfaceMaterialFactory::create(
        *axis, std::forward<accessor_t>(accessor),
        GridAccessJsonConverter::boundToGridLocal1DimDelegateFromJson(
            jMaterial.at("bound_to_grid_local")),
        GridAccessJsonConverter::globalToGridLocal1DimDelegateFromJson(
            jMaterial.at("global_to_grid_local")),
        payload);
  }
  if (jAxes.size() == 2u) {
    std::unique_ptr<IAxis> axis0 = AxisJsonConverter::fromJson(jAxes.at(0u));
    std::unique_ptr<IAxis> axis1 = AxisJsonConverter::fromJson(jAxes.at(1u));
    auto payload = gridPayload2D<typename accessor_t::grid_value_type>(
        jGrid, axis0->getNBins(), axis1->getNBins(), fill);
    return GridSurfaceMaterialFactory::create(
        *axis0, *axis1, std::forward<accessor_t>(accessor),
        GridAccessJsonConverter::boundToGridLocal2DimDelegateFromJson(
            jMaterial.at("bound_to_grid_local")),
        GridAccessJsonConverter::globalToGridLocal2DimDelegateFromJson(
            jMaterial.at("global_to_grid_local")),
        payload);
  }
  throw std::invalid_argument(
      "SurfaceMaterialJsonConverter: only 1D and 2D material grids are "
      "supported");
}

std::vector<MaterialSlab> slabsFromJson(const nlohmann::json& jSlabs) {
  std::vector<MaterialSlab> slabs;
  slabs.reserve(jSlabs.size());
  for (const auto& jSlab : jSlabs) {
    MaterialSlab slab = MaterialSlab::Nothing();
    from_json(jSlab, slab);
    slabs.push_back(slab);
  }
  return slabs;
}

std::unique_ptr<const ISurfaceMaterial> gridFromJson(
    const nlohmann::json& jMaterial, const DecodeContext& ctx) {
  const nlohmann::json& jAccessor = jMaterial.at("accessor");
  const nlohmann::json& jGrid = jAccessor.at("grid");
  std::string accessorType = jAccessor.at("type").get<std::string>();

  if (accessorType == kIndexedAccessorTag) {
    return createGridMaterial(
        jMaterial, jGrid,
        IndexedMaterialAccessor{slabsFromJson(jAccessor.at("storage_vector"))},
        std::size_t{0u});
  }
  if (accessorType == kGloballyIndexedAccessorTag) {
    bool sharedEntries = jAccessor.value("shared_entries", false);
    MaterialSlabStore store;
    if (jAccessor.contains("store")) {
      // Resolved through the document store table, so that grids referencing
      // the same id keep sharing one allocation
      store = ctx.store(jAccessor.at("store").get<std::size_t>());
    } else {
      store = std::make_shared<std::vector<MaterialSlab>>(
          slabsFromJson(jAccessor.at("storage_vector")));
    }
    return createGridMaterial(
        jMaterial, jGrid,
        GloballyIndexedMaterialAccessor{std::move(store), sharedEntries},
        std::size_t{0u});
  }
  if (accessorType == kGridAccessorTag) {
    return createGridMaterial(jMaterial, jGrid, GridMaterialAccessor{},
                              MaterialSlab::Nothing());
  }
  throw std::invalid_argument(
      "SurfaceMaterialJsonConverter: unsupported grid material accessor: " +
      accessorType);
}

}  // namespace

Acts::SurfaceMaterialJsonConverter::Config
Acts::SurfaceMaterialJsonConverter::Config::defaultConfig() {
  Config cfg;

  cfg.encoder.registerFunction(homogeneousToJson);
  cfg.encoder.registerFunction(binnedToJson);
  cfg.encoder.registerFunction(protoToJson);
  cfg.encoder.registerFunction(protoGridToJson);
  cfg.encoder.registerFunction(mergedMarkerToJson);
  // The whole grid material family is covered by the two abstract bases,
  // which keeps the concrete grid shapes out of this translation unit
  cfg.encoder.registerFunction(slabGridMaterialToJson);
  cfg.encoder.registerFunction(indexGridMaterialToJson);

  cfg.decoder.registerKind(kHomogeneousTag, homogeneousFromJson);
  cfg.decoder.registerKind(kBinnedTag, binnedFromJson);
  cfg.decoder.registerKind(kProtoTag, protoFromJson);
  cfg.decoder.registerKind(kProtoGridTag, protoGridFromJson);
  cfg.decoder.registerKind(kMergedMarkerTag, mergedMarkerFromJson);
  cfg.decoder.registerKind(kGridTag, gridFromJson);

  return cfg;
}

const Acts::SurfaceMaterialJsonConverter::Config&
Acts::SurfaceMaterialJsonConverter::defaultConfig() {
  static const Config cfg = Config::defaultConfig();
  return cfg;
}

nlohmann::json Acts::SurfaceMaterialJsonConverter::toJson(
    const ISurfaceMaterial& material, EncodeContext& context,
    const Config& config) {
  return config.encoder(material, context);
}

nlohmann::json Acts::SurfaceMaterialJsonConverter::toJson(
    const ISurfaceMaterial& material, const Config& config) {
  // A context without a store table makes the encoders inline their stores
  EncodeContext context;
  return toJson(material, context, config);
}

std::unique_ptr<const Acts::ISurfaceMaterial>
Acts::SurfaceMaterialJsonConverter::fromJson(const nlohmann::json& jMaterial,
                                             const DecodeContext& context,
                                             const Config& config) {
  // Surfaces that are flagged out of the mapping carry no material
  if (jMaterial.contains(jsonKey().mapkey) &&
      jMaterial.at(jsonKey().mapkey) == false) {
    return nullptr;
  }
  return config.decoder(jMaterial, context);
}

std::unique_ptr<const Acts::ISurfaceMaterial>
Acts::SurfaceMaterialJsonConverter::fromJson(const nlohmann::json& jMaterial,
                                             const Config& config) {
  const DecodeContext context;
  return fromJson(jMaterial, context, config);
}

void Acts::to_json(nlohmann::json& j,
                   const std::shared_ptr<const ISurfaceMaterial>& material) {
  if (material == nullptr) {
    return;
  }
  j[jsonKey().materialkey] = SurfaceMaterialJsonConverter::toJson(*material);
}

void Acts::from_json(const nlohmann::json& j,
                     std::shared_ptr<const ISurfaceMaterial>& material) {
  material = nullptr;
  if (!j.contains(jsonKey().materialkey) ||
      j.at(jsonKey().materialkey).is_null()) {
    return;
  }
  material =
      SurfaceMaterialJsonConverter::fromJson(j.at(jsonKey().materialkey));
}

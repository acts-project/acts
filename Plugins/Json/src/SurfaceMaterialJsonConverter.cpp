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
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/MergedMaterialMarker.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"
#include "ActsPlugins/Json/AxisSpecJsonConverter.hpp"
#include "ActsPlugins/Json/GridJsonConverter.hpp"
#include "ActsPlugins/Json/MaterialJsonConverter.hpp"
#include "ActsPlugins/Json/UtilitiesJsonConverter.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
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
constexpr const char* kDirectAccessorTag = "direct";
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

/// Convert the axes of a resolved multi-axis
nlohmann::json axesToJson(const IMultiAxis2D& multiAxis) {
  nlohmann::json jAxes = nlohmann::json::array();
  for (std::size_t ia = 0u; ia < multiAxis.getNAxes(); ++ia) {
    jAxes.push_back(AxisJsonConverter::toJson(multiAxis.getAxis(ia)));
  }
  return jAxes;
}

/// Write the per bin payload as [local bins, value] entries over the regular
/// bins. The local bins are 1-based, the storage is addressed by global bin.
template <typename value_at_t>
nlohmann::json gridDataToJson(const IMultiAxis2D& multiAxis,
                              value_at_t&& valueAt) {
  nlohmann::json jData = nlohmann::json::array();
  IMultiAxis2D::LocalBins nBins = multiAxis.getNBins();
  for (std::size_t ib0 = 1u; ib0 <= nBins[0u]; ++ib0) {
    for (std::size_t ib1 = 1u; ib1 <= nBins[1u]; ++ib1) {
      IMultiAxis2D::LocalBins lBin{ib0, ib1};
      jData.push_back(nlohmann::json::array(
          {std::array<std::size_t, 2u>{ib0, ib1},
           valueAt(multiAxis.getGlobalBinFromLocalBins(lBin))}));
    }
  }
  return jData;
}

nlohmann::json slabsToJson(const std::vector<MaterialSlab>& slabs) {
  nlohmann::json jSlabs = nlohmann::json::array();
  for (const auto& msl : slabs) {
    jSlabs.push_back(nlohmann::json(msl));
  }
  return jSlabs;
}

/// Encoder for the whole grid material family. The storage backend is
/// resolved at runtime through the storage variant, so the concrete axis
/// types never appear here.
nlohmann::json gridMaterialToJson(const GridSurfaceMaterial& material,
                                  EncodeContext& ctx) {
  nlohmann::json jMaterial;
  jMaterial[jsonKey().typekey] = kGridTag;
  jMaterial[jsonKey().mapkey] = true;

  const IMultiAxis2D& multiAxis = material.multiAxis();
  nlohmann::json jGrid;
  jGrid["axes"] = axesToJson(multiAxis);

  nlohmann::json jAccessor;
  std::visit(
      [&]<typename storage_t>(const storage_t& storage) {
        if constexpr (std::is_same_v<storage_t, GridSurfaceMaterial::Direct>) {
          jAccessor["type"] = kDirectAccessorTag;
          jGrid["data"] = gridDataToJson(multiAxis, [&](std::size_t bin) {
            return nlohmann::json(storage.at(bin));
          });
        } else {
          constexpr bool isLocal =
              std::is_same_v<storage_t, GridSurfaceMaterial::Indexed>;
          jAccessor["type"] =
              isLocal ? kIndexedAccessorTag : kGloballyIndexedAccessorTag;
          // The slab store travels with the payload, so that a globally
          // indexed grid can be read back on its own
          if constexpr (isLocal) {
            jAccessor["storage_vector"] = slabsToJson(storage.material);
          } else {
            if (storage.material == nullptr) {
              throw std::invalid_argument(
                  "SurfaceMaterialJsonConverter: globally indexed material "
                  "without a slab store");
            }
            if (ctx.storeTableEnabled()) {
              // The store lives once in the document, the entry only
              // references it
              jAccessor["store"] = ctx.storeId(storage.material);
            } else {
              // Standalone payload, inline the store to keep it
              // self-contained
              jAccessor["storage_vector"] = slabsToJson(*storage.material);
            }
          }
          jGrid["data"] = gridDataToJson(multiAxis, [&](std::size_t bin) {
            return nlohmann::json(storage.indices.at(bin));
          });
        }
      },
      material.storage());

  jAccessor["grid"] = std::move(jGrid);
  jMaterial["accessor"] = std::move(jAccessor);
  return jMaterial;
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

/// Read the 2D grid payload in column major order, i.e. [i0][i1], from the
/// 1-based local bins the writer emits
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
  const nlohmann::json& jAxes = jGrid.at("axes");
  std::string accessorType = jAccessor.at("type").get<std::string>();

  if (jAxes.size() != 2u) {
    throw std::invalid_argument(
        "SurfaceMaterialJsonConverter: grid material needs exactly two axes");
  }
  std::unique_ptr<IAxis> axis0 = AxisJsonConverter::fromJson(jAxes.at(0u));
  std::unique_ptr<IAxis> axis1 = AxisJsonConverter::fromJson(jAxes.at(1u));

  if (accessorType == kDirectAccessorTag) {
    return GridSurfaceMaterial::createDirect(
        *axis0, *axis1,
        gridPayload2D<MaterialSlab>(jGrid, axis0->getNBins(), axis1->getNBins(),
                                    MaterialSlab::Nothing()));
  }

  auto indices = gridPayload2D<std::size_t>(jGrid, axis0->getNBins(),
                                            axis1->getNBins(), std::size_t{0u});
  if (accessorType == kIndexedAccessorTag) {
    return GridSurfaceMaterial::createIndexed(
        *axis0, *axis1, slabsFromJson(jAccessor.at("storage_vector")), indices);
  }
  if (accessorType == kGloballyIndexedAccessorTag) {
    MaterialSlabStore store;
    if (jAccessor.contains("store")) {
      // Resolved through the document store table, so that grids referencing
      // the same id keep sharing one allocation
      store = ctx.store(jAccessor.at("store").get<std::size_t>());
    } else {
      store = std::make_shared<std::vector<MaterialSlab>>(
          slabsFromJson(jAccessor.at("storage_vector")));
    }
    return GridSurfaceMaterial::createGloballyIndexed(
        *axis0, *axis1, std::move(store), indices);
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
  // One concrete class covers the whole grid material family, the storage
  // backend is a runtime variant rather than a template parameter
  cfg.encoder.registerFunction(gridMaterialToJson);

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
  // Without a document context the encoders inline their slab stores
  EncodeContext inlineContext;
  return toJson(material, inlineContext, config);
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
  // Without a document context a payload referencing a store is an error
  const DecodeContext emptyContext;
  return fromJson(jMaterial, emptyContext, config);
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

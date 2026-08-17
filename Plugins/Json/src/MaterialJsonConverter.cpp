// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/MaterialJsonConverter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterialFactory.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/InterpolatedMaterialMap.hpp"
#include "Acts/Material/MaterialGridHelper.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/MergedMaterialMarker.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "ActsPlugins/Json/GeometryJsonKeys.hpp"
#include "ActsPlugins/Json/GridJsonConverter.hpp"
#include "ActsPlugins/Json/UtilitiesJsonConverter.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

/// @brief Convert the axes of a resolved multi-axis to json
nlohmann::json axesToJson(const Acts::IMultiAxis2D& multiAxis) {
  nlohmann::json jAxes;
  for (std::size_t i = 0; i < multiAxis.getNAxes(); ++i) {
    jAxes.push_back(Acts::AxisJsonConverter::toJson(multiAxis.getAxis(i)));
  }
  return jAxes;
}

/// @brief Convert the per-bin payload of a resolved multi-axis to json,
/// shaped as [localBins, value] entries over the regular (non-overflow) bins
///
/// @tparam value_at_t callable, global bin index -> json-convertible value
template <typename value_at_t>
nlohmann::json gridDataToJson(const Acts::IMultiAxis2D& multiAxis,
                              value_at_t&& valueAt) {
  nlohmann::json jData;
  Acts::IMultiAxis2D::LocalBins nBins = multiAxis.getNBins();
  for (std::size_t i0 = 1; i0 <= nBins[0]; ++i0) {
    for (std::size_t i1 = 1; i1 <= nBins[1]; ++i1) {
      Acts::IMultiAxis2D::LocalBins lbin{i0, i1};
      std::size_t bin = multiAxis.getGlobalBinFromLocalBins(lbin);
      std::array<std::size_t, 2u> jLbin{i0, i1};
      jData.push_back(nlohmann::json::array({jLbin, valueAt(bin)}));
    }
  }
  return jData;
}

/// @brief Convert a @c GridSurfaceMaterial to json
///
/// Works generically for all 3 storage backends (direct / locally indexed /
/// globally indexed) via @c GridSurfaceMaterial::storage(), so it does not
/// need to enumerate concrete axis type combinations, nor dynamic_cast to a
/// per-backend wrapper type.
///
/// @param jMaterial the json object to write into
/// @param gridMaterial the grid surface material to convert
void writeGridSurfaceMaterial(nlohmann::json& jMaterial,
                              const Acts::GridSurfaceMaterial& gridMaterial) {
  jMaterial[Acts::jsonKey().typekey] = "grid";

  const Acts::IMultiAxis2D& multiAxis = gridMaterial.multiAxis();

  nlohmann::json jGrid;
  jGrid["axes"] = axesToJson(multiAxis);

  nlohmann::json jMaterialAccessor;
  std::visit(
      [&]<typename T>(const T& storage) {
        if constexpr (std::is_same_v<T, Acts::GridSurfaceMaterial::Direct>) {
          jMaterialAccessor["type"] = "direct";
          jGrid["data"] = gridDataToJson(multiAxis, [&](std::size_t bin) {
            return nlohmann::json(storage.at(bin));
          });
        } else if constexpr (std::is_same_v<
                                 T, Acts::GridSurfaceMaterial::Indexed>) {
          jMaterialAccessor["type"] = "indexed";
          nlohmann::json jMaterialData;
          for (const auto& msl : storage.material) {
            jMaterialData.push_back(msl);
          }
          jMaterialAccessor["storage_vector"] = jMaterialData;
          jGrid["data"] = gridDataToJson(multiAxis, [&](std::size_t bin) {
            return nlohmann::json(storage.indices.at(bin));
          });
        } else {
          jMaterialAccessor["type"] = "globally_indexed";
          jGrid["data"] = gridDataToJson(multiAxis, [&](std::size_t bin) {
            return nlohmann::json(storage.indices.at(bin));
          });
        }
      },
      gridMaterial.storage());

  jMaterialAccessor["grid"] = jGrid;
  jMaterial["accessor"] = jMaterialAccessor;
}

/// @brief Reconstruct a 2D grid payload from the json "data" array
///
/// @tparam value_type the grid payload type (MaterialSlab or std::size_t)
/// @param jData the json "data" array, entries shaped as [localBins, value]
/// @param nBins0 the number of bins along axis 0
/// @param nBins1 the number of bins along axis 1
template <typename value_type>
std::vector<std::vector<value_type>> readGridPayload2D(
    const nlohmann::json& jData, std::size_t nBins0, std::size_t nBins1) {
  std::vector<std::vector<value_type>> payload(nBins0,
                                               std::vector<value_type>(nBins1));
  for (const auto& jd : jData) {
    std::array<std::size_t, 2u> lbin = jd[0u];
    if (!jd[1u].is_null()) {
      payload[lbin[0u] - 1u][lbin[1u] - 1u] = jd[1u].get<value_type>();
    }
  }
  return payload;
}

/// @brief Read the locally indexed material vector from the json accessor
///
/// @param jMaterialAccessor the json "accessor" object
/// @return the material vector, in storage order
std::vector<Acts::MaterialSlab> readStorageVector(
    const nlohmann::json& jMaterialAccessor) {
  std::vector<Acts::MaterialSlab> materialVector;
  for (const auto& msl : jMaterialAccessor["storage_vector"]) {
    Acts::MaterialSlab mat = Acts::MaterialSlab::Nothing();
    from_json(msl, mat);
    materialVector.push_back(mat);
  }
  return materialVector;
}

/// @brief Reconstruct a @c GridSurfaceMaterial from json
///
/// Works generically for arbitrary axis type combinations via the
/// @c IAxis-based @c GridSurfaceMaterialFactory, so it is not limited to a
/// hardcoded set of equidistant bound/closed axis combinations. The grid is
/// always 2D, matching @c GridSurfaceMaterial's requirement.
///
/// @param jMaterial the json object to read from
/// @return a newly allocated surface material, or nullptr if unsupported
/// @throws std::runtime_error for "globally_indexed" material: the shared
///         material vector has no source at single-surface json scope
Acts::ISurfaceMaterial* gridSurfaceMaterialFromJson(nlohmann::json& jMaterial) {
  nlohmann::json jMaterialAccessor = jMaterial["accessor"];
  std::string accessorType = jMaterialAccessor["type"];

  if (accessorType == "globally_indexed") {
    throw std::runtime_error(
        "MaterialJsonConverter: reading a globally indexed "
        "GridSurfaceMaterial from json is not supported - the shared "
        "material vector has no source at single-surface json scope.");
  }

  nlohmann::json jGrid = jMaterialAccessor["grid"];
  nlohmann::json jGridAxes = jGrid["axes"];
  nlohmann::json jData = jGrid["data"];

  std::vector<std::unique_ptr<Acts::IAxis>> axes;
  for (const auto& jAxis : jGridAxes) {
    axes.push_back(Acts::AxisJsonConverter::fromJson(jAxis));
  }

  if (axes.size() != 2u) {
    return nullptr;
  }

  if (accessorType == "direct") {
    auto payload = readGridPayload2D<Acts::MaterialSlab>(
        jData, axes[0]->getNBins(), axes[1]->getNBins());
    return Acts::GridSurfaceMaterialFactory::createDirect(*axes[0], *axes[1],
                                                          payload)
        .release();
  }
  if (accessorType == "indexed") {
    auto payload = readGridPayload2D<std::size_t>(jData, axes[0]->getNBins(),
                                                  axes[1]->getNBins());
    return Acts::GridSurfaceMaterialFactory::createIndexed(
               *axes[0], *axes[1], readStorageVector(jMaterialAccessor),
               payload)
        .release();
  }

  return nullptr;
}

}  // namespace

void Acts::to_json(nlohmann::json& j, const Material& t) {
  if (t.isVacuum()) {
    return;
  }
  for (unsigned i = 0; i < t.parameters().size(); ++i) {
    j.push_back(t.parameters()[i]);
  }
}

void Acts::from_json(const nlohmann::json& j, Material& t) {
  if (j.is_null()) {
    return;
  }
  Acts::Material::ParametersVector params =
      Acts::Material::ParametersVector::Zero();
  for (auto i = params.size(); 0 < i--;) {
    // .at(...) ensures bound checks
    params[i] = j.at(i);
  }
  t = Acts::Material(params);
  return;
}

void Acts::to_json(nlohmann::json& j, const MaterialSlab& t) {
  nlohmann::json jmat(t.material());
  j["material"] = jmat;
  j["thickness"] = t.thickness();
}

void Acts::from_json(const nlohmann::json& j, MaterialSlab& t) {
  Material mat = Material::Vacuum();
  from_json(j.at("material"), mat);
  t = Acts::MaterialSlab(mat, j.at("thickness").get<float>());
}

void Acts::from_json(const nlohmann::json& j, MaterialSlabMatrix& t) {
  // the input data must be array[array[object]]
  for (auto& outer : j) {
    Acts::MaterialSlabVector mpVector;
    for (auto& inner : outer) {
      MaterialSlab mat = MaterialSlab::Nothing();
      from_json(inner, mat);
      mpVector.emplace_back(mat);
    }
    t.push_back(std::move(mpVector));
  }
}

void Acts::to_json(nlohmann::json& j, const surfaceMaterialPointer& material) {
  nlohmann::json jMaterial;
  // A bin utility needs to be written
  const Acts::BinUtility* bUtility = nullptr;

  // Marker material left behind by a lossy portal merge. It carries no actual
  // material, so only the type tag is written.
  if (dynamic_cast<const Acts::MergedMaterialMarker*>(material) != nullptr) {
    jMaterial[Acts::jsonKey().typekey] = "merged-material-marker";
    // Flag as "mapped" so the reader does not discard it.
    jMaterial[Acts::jsonKey().mapkey] = true;
    j[Acts::jsonKey().materialkey] = jMaterial;
    return;
  }

  // First: Check if we have a proto material
  auto psMaterial = dynamic_cast<const Acts::ProtoSurfaceMaterial*>(material);
  if (psMaterial != nullptr) {
    // Type is proto material
    jMaterial[Acts::jsonKey().typekey] = "proto";
    // Set mapping type
    nlohmann::json mapType(material->mappingType());
    jMaterial[Acts::jsonKey().maptype] = mapType;
    // by default the protoMaterial is not used for mapping
    jMaterial[Acts::jsonKey().mapkey] = false;
    // write the bin utility
    bUtility = &(psMaterial->binning());
    // Check in the number of bin is different from 1
    auto& binningData = bUtility->binningData();
    for (std::size_t ibin = 0; ibin < binningData.size(); ++ibin) {
      if (binningData[ibin].bins() > 1) {
        jMaterial[Acts::jsonKey().mapkey] = true;
        break;
      }
    }
    nlohmann::json jBin(*bUtility);
    jMaterial[Acts::jsonKey().binkey] = jBin;
    j[Acts::jsonKey().materialkey] = jMaterial;
    return;
  }

  // Second: check if we have a homogeneous material
  auto hsMaterial =
      dynamic_cast<const Acts::HomogeneousSurfaceMaterial*>(material);
  if (hsMaterial != nullptr) {
    // type is homogeneous
    jMaterial[Acts::jsonKey().typekey] = "homogeneous";
    // Set mapping type
    nlohmann::json mapType(material->mappingType());
    jMaterial[Acts::jsonKey().maptype] = mapType;
    // Material has been mapped
    jMaterial[Acts::jsonKey().mapkey] = true;
    nlohmann::json jmat(hsMaterial->materialSlab());
    jMaterial[Acts::jsonKey().datakey] = nlohmann::json::array({
        nlohmann::json::array({
            jmat,
        }),
    });
    j[Acts::jsonKey().materialkey] = jMaterial;
    return;
  }

  // Next option remaining: BinnedSurface material
  auto bsMaterial = dynamic_cast<const Acts::BinnedSurfaceMaterial*>(material);
  if (bsMaterial != nullptr) {
    // type is binned
    jMaterial[Acts::jsonKey().typekey] = "binned";
    // Set mapping type
    nlohmann::json mapType(material->mappingType());
    jMaterial[Acts::jsonKey().maptype] = mapType;
    // Material has been mapped
    jMaterial[Acts::jsonKey().mapkey] = true;
    bUtility = &(bsMaterial->binUtility());
    // convert the data
    // get the material matrix
    nlohmann::json mmat = nlohmann::json::array();
    for (const auto& mpVector : bsMaterial->fullMaterial()) {
      nlohmann::json mvec = nlohmann::json::array();
      for (const auto& mp : mpVector) {
        nlohmann::json jmat(mp);
        mvec.push_back(jmat);
      }
      mmat.push_back(std::move(mvec));
    }
    jMaterial[Acts::jsonKey().datakey] = std::move(mmat);
    // write the bin utility
    nlohmann::json jBin(*bUtility);
    jMaterial[Acts::jsonKey().binkey] = jBin;
    j[Acts::jsonKey().materialkey] = jMaterial;
    return;
  }

  // Grid-based surface material: direct, indexed and globally indexed
  // storage all share the same I/O, branching only on the storage variant
  if (auto gridMaterial =
          dynamic_cast<const Acts::GridSurfaceMaterial*>(material);
      gridMaterial != nullptr) {
    writeGridSurfaceMaterial(jMaterial, *gridMaterial);
    j[Acts::jsonKey().materialkey] = jMaterial;
    return;
  }

  // No material the json object is left empty.
  return;
}

void Acts::from_json(const nlohmann::json& j,
                     surfaceMaterialPointer& material) {
  if (j.find(Acts::jsonKey().materialkey) == j.end()) {
    return;
  }
  nlohmann::json jMaterial = j[Acts::jsonKey().materialkey];
  // By default no material is return.
  material = nullptr;
  if (jMaterial[Acts::jsonKey().mapkey] == false) {
    return;
  }

  // Marker material left behind by a lossy portal merge
  if (jMaterial.contains(Acts::jsonKey().typekey) &&
      jMaterial[Acts::jsonKey().typekey] == "merged-material-marker") {
    material = new Acts::MergedMaterialMarker();
    return;
  }

  // Grid based material maps
  if (jMaterial[Acts::jsonKey().typekey] == "grid") {
    material = gridSurfaceMaterialFromJson(jMaterial);
    return;
  }

  // The bin utility and material
  Acts::BinUtility bUtility;
  Acts::MaterialSlabMatrix mpMatrix;
  Acts::MappingType mapType = Acts::MappingType::Default;
  for (auto& [key, value] : jMaterial.items()) {
    if (key == Acts::jsonKey().binkey && !value.empty()) {
      from_json(value, bUtility);
    }
    if (key == Acts::jsonKey().datakey && !value.empty()) {
      from_json(value, mpMatrix);
    }
    if (key == Acts::jsonKey().maptype && !value.empty()) {
      from_json(value, mapType);
    }
  }
  // Return the appropriate typr of material
  if (mpMatrix.empty()) {
    material = new Acts::ProtoSurfaceMaterial(bUtility, mapType);
  } else if (bUtility.bins() == 1) {
    material = new Acts::HomogeneousSurfaceMaterial(mpMatrix[0][0], 1, mapType);
  } else {
    material = new Acts::BinnedSurfaceMaterial(bUtility, mpMatrix, 1, mapType);
  }
}

void Acts::to_json(nlohmann::json& j, const volumeMaterialPointer& material) {
  nlohmann::json jMaterial;
  // A bin utility needs to be written
  const Acts::BinUtility* bUtility = nullptr;
  // Check if we have a proto material
  auto pvMaterial = dynamic_cast<const Acts::ProtoVolumeMaterial*>(material);
  if (pvMaterial != nullptr) {
    // Type is proto material
    jMaterial[Acts::jsonKey().typekey] = "proto";
    // By default the protoMaterial is not used for mapping
    jMaterial[Acts::jsonKey().mapkey] = false;
    bUtility = &(pvMaterial->binUtility());
    // Check in the number of bin is different from 1
    auto& binningData = bUtility->binningData();
    for (std::size_t ibin = 0; ibin < binningData.size(); ++ibin) {
      if (binningData[ibin].bins() > 1) {
        jMaterial[Acts::jsonKey().mapkey] = true;
        break;
      }
    }
    // Write the bin utility
    nlohmann::json jBin(*bUtility);
    jMaterial[Acts::jsonKey().binkey] = jBin;
    j[Acts::jsonKey().materialkey] = jMaterial;
    return;
  }
  // Now check if we have a homogeneous material
  auto hvMaterial =
      dynamic_cast<const Acts::HomogeneousVolumeMaterial*>(material);
  if (hvMaterial != nullptr) {
    // type is homogeneous
    jMaterial[Acts::jsonKey().typekey] = "homogeneous";
    jMaterial[Acts::jsonKey().mapkey] = true;
    // array of encoded materials w/ one entry
    nlohmann::json jmat(hvMaterial->material({0, 0, 0}));
    jMaterial[Acts::jsonKey().datakey] = nlohmann::json::array({
        jmat,
    });
    j[Acts::jsonKey().materialkey] = jMaterial;
    return;
  }
  // Only option remaining: material map
  auto bvMaterial2D = dynamic_cast<const Acts::InterpolatedMaterialMap<
      Acts::MaterialMapLookup<Acts::MaterialGrid2D>>*>(material);
  // Now check if we have a 2D map
  if (bvMaterial2D != nullptr) {
    // type is binned
    jMaterial[Acts::jsonKey().typekey] = "interpolated2D";
    jMaterial[Acts::jsonKey().mapkey] = true;
    bUtility = &(bvMaterial2D->binUtility());
    // convert the data
    nlohmann::json mmat = nlohmann::json::array();
    Acts::MaterialGrid2D grid = bvMaterial2D->getMapper().getGrid();
    for (std::size_t bin = 0; bin < grid.size(); bin++) {
      nlohmann::json jmat(Material(grid.at(bin)));
      mmat.push_back(jmat);
    }
    jMaterial[Acts::jsonKey().datakey] = std::move(mmat);
    // Write the bin utility
    nlohmann::json jBin(*bUtility);
    jMaterial[Acts::jsonKey().binkey] = jBin;
    j[Acts::jsonKey().materialkey] = jMaterial;
    return;
  }
  // Only option remaining: material map
  auto bvMaterial3D = dynamic_cast<const Acts::InterpolatedMaterialMap<
      Acts::MaterialMapLookup<Acts::MaterialGrid3D>>*>(material);
  // Now check if we have a 3D map
  if (bvMaterial3D != nullptr) {
    // type is binned
    jMaterial[Acts::jsonKey().typekey] = "interpolated3D";
    jMaterial[Acts::jsonKey().mapkey] = true;
    bUtility = &(bvMaterial3D->binUtility());
    // convert the data
    nlohmann::json mmat = nlohmann::json::array();
    Acts::MaterialGrid3D grid = bvMaterial3D->getMapper().getGrid();
    for (std::size_t bin = 0; bin < grid.size(); bin++) {
      nlohmann::json jmat(Material(grid.at(bin)));
      mmat.push_back(jmat);
    }
    jMaterial[Acts::jsonKey().datakey] = std::move(mmat);
    // Write the bin utility
    nlohmann::json jBin(*bUtility);
    jMaterial[Acts::jsonKey().binkey] = jBin;
    j[Acts::jsonKey().materialkey] = jMaterial;
    return;
  }
}

void Acts::from_json(const nlohmann::json& j, volumeMaterialPointer& material) {
  if (j.find(Acts::jsonKey().materialkey) == j.end()) {
    return;
  }
  nlohmann::json jMaterial = j[Acts::jsonKey().materialkey];
  // By default no material is return.
  material = nullptr;
  if (jMaterial[Acts::jsonKey().mapkey] == false) {
    return;
  }
  // The bin utility and material
  Acts::BinUtility bUtility;
  std::vector<Acts::Material> mmat;
  for (auto& [key, value] : jMaterial.items()) {
    if (key == Acts::jsonKey().binkey && !value.empty()) {
      from_json(value, bUtility);
    }
    if (key == Acts::jsonKey().datakey && !value.empty()) {
      for (const auto& bin : value) {
        Acts::Material mat = Material::Vacuum();
        from_json(bin, mat);
        mmat.push_back(mat);
      }
    }
  }
  // We have protoMaterial
  if (mmat.empty()) {
    material = new Acts::ProtoVolumeMaterial(bUtility);
    return;
  }
  if (mmat.size() == 1) {
    material = new Acts::HomogeneousVolumeMaterial(mmat[0]);
    return;
  }
  if (bUtility.dimensions() == 2) {
    std::function<Acts::Vector2(Acts::Vector3)> transfoGlobalToLocal;
    Acts::Grid2D grid = createGrid2D(bUtility, transfoGlobalToLocal);

    Acts::Grid2D::point_t min = grid.multiAxis().getMinPoint();
    Acts::Grid2D::point_t max = grid.multiAxis().getMaxPoint();
    Acts::Grid2D::index_t nBins = grid.multiAxis().getNBins();

    Acts::EAxis axis1(min[0], max[0], nBins[0]);
    Acts::EAxis axis2(min[1], max[1], nBins[1]);

    // Build the grid and fill it with data
    Acts::MaterialGrid2D mGrid(std::make_tuple(axis1, axis2));

    for (std::size_t bin = 0; bin < mmat.size(); bin++) {
      mGrid.at(bin) = mmat[bin].parameters();
    }
    Acts::MaterialMapLookup<Acts::MaterialGrid2D> matMap(transfoGlobalToLocal,
                                                         mGrid);
    material = new Acts::InterpolatedMaterialMap<
        Acts::MaterialMapLookup<Acts::MaterialGrid2D>>(std::move(matMap),
                                                       bUtility);
    return;
  }
  if (bUtility.dimensions() == 3) {
    std::function<Acts::Vector3(Acts::Vector3)> transfoGlobalToLocal;
    Acts::Grid3D grid = createGrid3D(bUtility, transfoGlobalToLocal);

    Acts::Grid3D::point_t min = grid.multiAxis().getMinPoint();
    Acts::Grid3D::point_t max = grid.multiAxis().getMaxPoint();
    Acts::Grid3D::index_t nBins = grid.multiAxis().getNBins();

    Acts::EAxis axis1(min[0], max[0], nBins[0]);
    Acts::EAxis axis2(min[1], max[1], nBins[1]);
    Acts::EAxis axis3(min[2], max[2], nBins[2]);

    // Build the grid and fill it with data
    Acts::MaterialGrid3D mGrid(std::make_tuple(axis1, axis2, axis3));

    for (std::size_t bin = 0; bin < mmat.size(); bin++) {
      mGrid.at(bin) = mmat[bin].parameters();
    }
    Acts::MaterialMapLookup<Acts::MaterialGrid3D> matMap(transfoGlobalToLocal,
                                                         mGrid);
    material = new Acts::InterpolatedMaterialMap<
        Acts::MaterialMapLookup<Acts::MaterialGrid3D>>(std::move(matMap),
                                                       bUtility);
    return;
  }
}

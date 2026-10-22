// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/MaterialJsonConverter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/GridVolumeMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "ActsPlugins/Json/AlgebraJsonConverter.hpp"
#include "ActsPlugins/Json/AxisSpecJsonConverter.hpp"
#include "ActsPlugins/Json/GeometryJsonKeys.hpp"
#include "ActsPlugins/Json/SurfaceMaterialJsonConverter.hpp"
#include "ActsPlugins/Json/UtilitiesJsonConverter.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

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
  if (material == nullptr) {
    return;
  }
  j[Acts::jsonKey().materialkey] =
      SurfaceMaterialJsonConverter::toJson(*material);
}

void Acts::from_json(const nlohmann::json& j,
                     surfaceMaterialPointer& material) {
  material = nullptr;
  if (j.find(Acts::jsonKey().materialkey) == j.end()) {
    return;
  }
  material =
      SurfaceMaterialJsonConverter::fromJson(j[Acts::jsonKey().materialkey])
          .release();
}

void Acts::to_json(nlohmann::json& j, const volumeMaterialPointer& material) {
  if (material == nullptr) {
    return;
  }
  nlohmann::json jMaterial;
  if (const auto* grid = dynamic_cast<const GridVolumeMaterial*>(material)) {
    jMaterial = {
        {"type", "grid"},
        {"mapMaterial", true},
        {"binning", MultiAxisSpecJsonConverter::toJson(grid->binning())},
        {"transform", Transform3JsonConverter::toJson(grid->transform())}};
    jMaterial["storage"] = std::visit(
        [](const auto& values) -> nlohmann::json {
          using T = std::decay_t<decltype(values)>;
          if constexpr (std::is_same_v<T, GridVolumeMaterial::Direct>) {
            return {{"type", "direct"}, {"material", values}};
          } else if constexpr (std::is_same_v<T, GridVolumeMaterial::Indexed>) {
            return {{"type", "indexed"},
                    {"indices", values.indices},
                    {"material", values.material}};
          } else {
            // Inline shared palettes to keep a volume payload self-contained.
            return {{"type", "globally-indexed"},
                    {"indices", values.indices},
                    {"material", *values.material}};
          }
        },
        grid->storage());
    j[Acts::jsonKey().materialkey] = std::move(jMaterial);
    return;
  }
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
}

void Acts::from_json(const nlohmann::json& j, volumeMaterialPointer& material) {
  if (j.find(Acts::jsonKey().materialkey) == j.end()) {
    return;
  }
  nlohmann::json jMaterial = j[Acts::jsonKey().materialkey];
  const auto type = jMaterial.value(Acts::jsonKey().typekey, std::string{});
  if (type == "interpolated2D" || type == "interpolated3D") {
    throw std::invalid_argument(
        "Interpolated volume material maps are no longer supported");
  }
  // By default no material is return.
  material = nullptr;
  if (jMaterial[Acts::jsonKey().mapkey] == false) {
    return;
  }
  if (type == "grid") {
    const auto spec =
        MultiAxisSpecJsonConverter::fromJson(jMaterial.at("binning"));
    if (spec.size() != 3) {
      throw std::invalid_argument("Grid volume material requires three axes");
    }
    MultiAxisSpec3D binning(
        std::array{spec.axisSpec(0), spec.axisSpec(1), spec.axisSpec(2)});
    Transform3 transform = Transform3::Identity();
    if (jMaterial.contains("transform")) {
      transform = Transform3JsonConverter::fromJson(jMaterial.at("transform"));
    }
    const auto& store = jMaterial.at("storage");
    const auto storageType = store.at("type").get<std::string>();
    std::vector<Material> palette;
    for (const auto& entry : store.at("material")) {
      auto value = Material::Vacuum();
      from_json(entry, value);
      palette.push_back(value);
    }
    GridVolumeMaterial::Storage storage;
    if (storageType == "direct") {
      storage = std::move(palette);
    } else if (storageType == "indexed") {
      storage = GridVolumeMaterial::Indexed{
          store.at("indices").get<std::vector<std::size_t>>(),
          std::move(palette)};
    } else if (storageType == "globally-indexed") {
      storage = GridVolumeMaterial::GloballyIndexed{
          store.at("indices").get<std::vector<std::size_t>>(),
          std::make_shared<std::vector<Material>>(std::move(palette))};
    } else {
      throw std::invalid_argument("Unknown grid volume material storage type");
    }
    material = new GridVolumeMaterial(std::move(binning), std::move(storage),
                                      transform);
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
  throw std::invalid_argument(
      "Interpolated volume material maps are no longer supported");
}

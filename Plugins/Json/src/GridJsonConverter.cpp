// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/GridJsonConverter.hpp"

#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Utilities/IAxis.hpp"

nlohmann::json Acts::AxisJsonConverter::toJson(const IAxis& ia) {
  nlohmann::json jAxis;

  jAxis["boundary_type"] = ia.getBoundaryType();
  // type, range, bins or boundaries
  if (ia.isEquidistant()) {
    jAxis["type"] = detail::AxisType::Equidistant;
    jAxis["range"] = std::array<ActsScalar, 2u>({ia.getMin(), ia.getMax()});
    jAxis["bins"] = ia.getNBins();
  } else {
    jAxis["type"] = detail::AxisType::Variable;
    jAxis["boundaries"] = ia.getBinEdges();
  }
  return jAxis;
}

nlohmann::json Acts::AxisJsonConverter::toJsonDetray(const IAxis& ia) {
  nlohmann::json jAxis;
  jAxis["bounds"] =
      ia.getBoundaryType() == Acts::detail::AxisBoundaryType::Bound ? 1 : 2;
  jAxis["binning"] = ia.isEquidistant() ? 0 : 1;
  jAxis["bins"] = ia.getNBins();
  if (ia.isEquidistant()) {
    std::array<ActsScalar, 2u> range = {ia.getBinEdges().front(),
                                        ia.getBinEdges().back()};
    jAxis["edges"] = range;

  } else {
    jAxis["edges"] = ia.getBinEdges();
  }
  return jAxis;
}

namespace {

template <typename Subspace>
void encodeSubspace(nlohmann::json& j,
                    const Acts::GridAccess::IGlobalToGridLocal& ga,
                    const Subspace& /*unused*/) {
  const Subspace* subspace = dynamic_cast<const Subspace*>(&ga);
  if (subspace != nullptr) {
    j["type"] = "subspace";
    j["accessors"] = subspace->bValues;
  }
}

template <typename Subspace>
void encodeTransformedSubspace(nlohmann::json& j,
                               const Acts::GridAccess::IGlobalToGridLocal& ga,
                               const Subspace& s) {
  const Acts::GridAccess::Affine3Transformed<Subspace>* tsubspace =
      dynamic_cast<const Acts::GridAccess::Affine3Transformed<Subspace>*>(&ga);
  if (tsubspace != nullptr) {
    encodeSubspace(j, tsubspace->globalToGridLocal, s);
    j["transform"] =
        Acts::Transform3JsonConverter::toJson(tsubspace->transform);
  }
}

template <typename... Args>
void encodeSubspaces(
    nlohmann::json& globalToGridLocalJson,
    const Acts::GridAccess::IGlobalToGridLocal& globalToGridLocal,
    bool transformed, const std::tuple<Args...>& tAcessors) {
  if (transformed) {
    std::apply(
        [&](auto&&... vals) {
          (encodeTransformedSubspace(globalToGridLocalJson, globalToGridLocal,
                                     vals),
           ...);
        },
        tAcessors);
  } else {
    std::apply(
        [&](auto&&... vals) {
          (encodeSubspace(globalToGridLocalJson, globalToGridLocal, vals), ...);
        },
        tAcessors);
  }
}

template <Acts::BinningValue... Args>
std::unique_ptr<Acts::GridAccess::IGlobalToGridLocal> decodeSubspace(
    const nlohmann::json& jga) {
  std::unique_ptr<Acts::GridAccess::IGlobalToGridLocal> ga = nullptr;
  if (jga.find("transform") != jga.end()) {
    Acts::Transform3 transform =
        Acts::Transform3JsonConverter::fromJson(jga.at("transform"));
    Acts::GridAccess::GlobalSubspace<Args...> s;
    ga = std::make_unique<Acts::GridAccess::Affine3Transformed<
        Acts::GridAccess::GlobalSubspace<Args...>>>(std::move(s), transform);
  } else {
    ga = std::make_unique<Acts::GridAccess::GlobalSubspace<Args...>>();
  }
  return ga;
}

}  // namespace

nlohmann::json Acts::GridAccessJsonConverter::toJson(
    const GridAccess::IGlobalToGridLocal& gridToGridLocal) {
  nlohmann::json gridToGridLocalJson;

  std::array<bool, 2u> transformOptions = {false, true};

  // One dimensional sub spaces
  const std::tuple<
      GridAccess::GlobalSubspace<binX>, GridAccess::GlobalSubspace<binY>,
      GridAccess::GlobalSubspace<binZ>, GridAccess::GlobalSubspace<binR>,
      GridAccess::GlobalSubspace<binPhi>, GridAccess::GlobalSubspace<binEta>>
      oneDimSubspaces = {};

  for (bool transform : transformOptions) {
    encodeSubspaces(gridToGridLocalJson, gridToGridLocal, transform,
                    oneDimSubspaces);
    if (!gridToGridLocalJson.empty()) {
      return gridToGridLocalJson;
    }
  }

  // Usfeul two dimensional sub spaces
  const std::tuple<GridAccess::GlobalSubspace<binX, binY>,
                   GridAccess::GlobalSubspace<binY, binX>,
                   GridAccess::GlobalSubspace<binX, binZ>,
                   GridAccess::GlobalSubspace<binZ, binX>,
                   GridAccess::GlobalSubspace<binY, binZ>,
                   GridAccess::GlobalSubspace<binZ, binY>,
                   GridAccess::GlobalSubspace<binR, binPhi>,
                   GridAccess::GlobalSubspace<binPhi, binR>,
                   GridAccess::GlobalSubspace<binZ, binPhi>,
                   GridAccess::GlobalSubspace<binPhi, binZ>>
      twoDimSubspaces = {};

  for (bool transform : transformOptions) {
    encodeSubspaces(gridToGridLocalJson, gridToGridLocal, transform,
                    twoDimSubspaces);
    if (!gridToGridLocalJson.empty()) {
      return gridToGridLocalJson;
    }
  }
  return gridToGridLocalJson;
}

std::unique_ptr<Acts::GridAccess::IGlobalToGridLocal>
Acts::GridAccessJsonConverter::globalToGridLocalFromJson(
    const nlohmann::json& gToGridLocalJson) {
  std::unique_ptr<Acts::GridAccess::IGlobalToGridLocal> gridToGridLocal =
      nullptr;

  std::vector<BinningValue> accessors =
      gToGridLocalJson.at("accessors").get<std::vector<BinningValue>>();

  // Switch and fill for 1D
  if (accessors.size() == 1u) {
    switch (accessors[0]) {
      case binX:
        gridToGridLocal = decodeSubspace<binX>(gToGridLocalJson);
        break;
      case binY:
        gridToGridLocal = decodeSubspace<binY>(gToGridLocalJson);
        break;
      case binZ:
        gridToGridLocal = decodeSubspace<binZ>(gToGridLocalJson);
        break;
      case binR:
        gridToGridLocal = decodeSubspace<binR>(gToGridLocalJson);
        break;
      case binPhi:
        gridToGridLocal = decodeSubspace<binPhi>(gToGridLocalJson);
        break;
      case binEta:
        gridToGridLocal = decodeSubspace<binEta>(gToGridLocalJson);
        break;
      default:
        break;
    }
  }

  // Switch and fill for 2D
  if (accessors.size() == 2u) {
    if (accessors == std::vector<BinningValue>{binX, binY}) {
      gridToGridLocal = decodeSubspace<binX, binY>(gToGridLocalJson);
    } else if (accessors == std::vector<BinningValue>{binY, binX}) {
      gridToGridLocal = decodeSubspace<binY, binX>(gToGridLocalJson);
    } else if (accessors == std::vector<BinningValue>{binX, binZ}) {
      gridToGridLocal = decodeSubspace<binX, binZ>(gToGridLocalJson);
    } else if (accessors == std::vector<BinningValue>{binZ, binX}) {
      gridToGridLocal = decodeSubspace<binZ, binX>(gToGridLocalJson);
    } else if (accessors == std::vector<BinningValue>{binY, binZ}) {
      gridToGridLocal = decodeSubspace<binY, binZ>(gToGridLocalJson);
    } else if (accessors == std::vector<BinningValue>{binZ, binY}) {
      gridToGridLocal = decodeSubspace<binZ, binY>(gToGridLocalJson);
    } else if (accessors == std::vector<BinningValue>{binR, binPhi}) {
      gridToGridLocal = decodeSubspace<binR, binPhi>(gToGridLocalJson);
    } else if (accessors == std::vector<BinningValue>{binPhi, binR}) {
      gridToGridLocal = decodeSubspace<binPhi, binR>(gToGridLocalJson);
    } else if (accessors == std::vector<BinningValue>{binZ, binPhi}) {
      gridToGridLocal = decodeSubspace<binZ, binPhi>(gToGridLocalJson);
    } else if (accessors == std::vector<BinningValue>{binPhi, binZ}) {
      gridToGridLocal = decodeSubspace<binPhi, binZ>(gToGridLocalJson);
    }
  }
  return gridToGridLocal;
}

nlohmann::json Acts::GridAccessJsonConverter::toJson(
    const GridAccess::IBoundToGridLocal& bToGridLocal) {
  nlohmann::json bToGridLocalJson;

  auto localSubSpace0 =
      dynamic_cast<const GridAccess::LocalSubspace<0u>*>(&bToGridLocal);
  if (localSubSpace0 != nullptr) {
    bToGridLocalJson["type"] = "subspace";
    bToGridLocalJson["accessors"] = localSubSpace0->accessors;
  }

  auto localSubSpace1 =
      dynamic_cast<const GridAccess::LocalSubspace<1u>*>(&bToGridLocal);
  if (localSubSpace1 != nullptr) {
    bToGridLocalJson["type"] = "subspace";
    bToGridLocalJson["accessors"] = localSubSpace1->accessors;
  }

  auto localSubSpace01 =
      dynamic_cast<const GridAccess::LocalSubspace<0u, 1u>*>(&bToGridLocal);
  if (localSubSpace01 != nullptr) {
    bToGridLocalJson["type"] = "subspace";
    bToGridLocalJson["accessors"] = localSubSpace01->accessors;
  }

  auto localSubSpace10 =
      dynamic_cast<const GridAccess::LocalSubspace<1u, 0u>*>(&bToGridLocal);
  if (localSubSpace10 != nullptr) {
    bToGridLocalJson["type"] = "subspace";
    bToGridLocalJson["accessors"] = localSubSpace10->accessors;
  }

  auto boundCylinderToZPhi =
      dynamic_cast<const GridAccess::BoundCylinderToZPhi*>(&bToGridLocal);
  if (boundCylinderToZPhi != nullptr) {
    bToGridLocalJson["type"] = "cylinder_to_zphi";
    bToGridLocalJson["radius"] = boundCylinderToZPhi->radius;
    bToGridLocalJson["shift"] = boundCylinderToZPhi->shift;
  }

  return bToGridLocalJson;
}

std::unique_ptr<Acts::GridAccess::IBoundToGridLocal>
Acts::GridAccessJsonConverter::boundToGridLocalFromJson(
    const nlohmann::json& bToGridLocalJson) {
  std::unique_ptr<Acts::GridAccess::IBoundToGridLocal> bToGridLocal = nullptr;
  std::string type = bToGridLocalJson.at("type").get<std::string>();
  if (type == "subspace") {
    std::vector<std::size_t> accessors =
        bToGridLocalJson.at("accessors").get<std::vector<std::size_t>>();
    if (accessors.size() == 1 && accessors[0] == 0) {
      bToGridLocal = std::make_unique<Acts::GridAccess::LocalSubspace<0u>>();
    } else if (accessors.size() == 1 && accessors[0] == 1) {
      bToGridLocal = std::make_unique<Acts::GridAccess::LocalSubspace<1u>>();
    } else if (accessors.size() == 2 && accessors[0] == 0 &&
               accessors[1] == 1) {
      bToGridLocal =
          std::make_unique<Acts::GridAccess::LocalSubspace<0u, 1u>>();
    } else if (accessors.size() == 2 && accessors[0] == 1 &&
               accessors[1] == 0) {
      bToGridLocal =
          std::make_unique<Acts::GridAccess::LocalSubspace<1u, 0u>>();
    }
  } else if (type == "cylinder_to_zphi") {
    ActsScalar radius = bToGridLocalJson.at("radius").get<ActsScalar>();
    ActsScalar shift = bToGridLocalJson.at("shift").get<ActsScalar>();
    bToGridLocal =
        std::make_unique<Acts::GridAccess::BoundCylinderToZPhi>(radius, shift);
  }
  return bToGridLocal;
}

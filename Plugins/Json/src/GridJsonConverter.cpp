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
void encodeSubspace(
    nlohmann::json& jGlobalToGridLocal,
    const Acts::GridAccess::IGlobalToGridLocal& globalToGridLocal,
    const Subspace& /*unused*/) {
  const Subspace* subspace = dynamic_cast<const Subspace*>(&globalToGridLocal);
  if (subspace != nullptr) {
    jGlobalToGridLocal["type"] = "subspace";
    jGlobalToGridLocal["accessors"] = subspace->bValues;
  }
}

template <typename Subspace>
void encodeTransformedSubspace(
    nlohmann::json& jGlobalToGridLocal,
    const Acts::GridAccess::IGlobalToGridLocal& globalToGridLocal,
    const Subspace& subscpace) {
  const Acts::GridAccess::Affine3Transformed<Subspace>* tsubspace =
      dynamic_cast<const Acts::GridAccess::Affine3Transformed<Subspace>*>(
          &globalToGridLocal);
  if (tsubspace != nullptr) {
    encodeSubspace(jGlobalToGridLocal, tsubspace->globalToGridLocal, subscpace);
    jGlobalToGridLocal["transform"] =
        Acts::Transform3JsonConverter::toJson(tsubspace->transform);
  }
}

template <typename... Args>
void encodeSubspaces(
    nlohmann::json& jGlobalToGridLocal,
    const Acts::GridAccess::IGlobalToGridLocal& globalToGridLocal,
    bool transformed, const std::tuple<Args...>& tAcessors) {
  if (transformed) {
    std::apply(
        [&](auto&&... vals) {
          (encodeTransformedSubspace(jGlobalToGridLocal, globalToGridLocal,
                                     vals),
           ...);
        },
        tAcessors);
  } else {
    std::apply(
        [&](auto&&... vals) {
          (encodeSubspace(jGlobalToGridLocal, globalToGridLocal, vals), ...);
        },
        tAcessors);
  }
}

template <Acts::BinningValue... Args>
std::unique_ptr<Acts::GridAccess::IGlobalToGridLocal> decodeSubspace(
    const nlohmann::json& jGlobalToGridLocal) {
  std::unique_ptr<Acts::GridAccess::IGlobalToGridLocal> globalToGridLocal =
      nullptr;
  if (jGlobalToGridLocal.find("transform") != jGlobalToGridLocal.end()) {
    Acts::Transform3 transform = Acts::Transform3JsonConverter::fromJson(
        jGlobalToGridLocal.at("transform"));
    Acts::GridAccess::GlobalSubspace<Args...> globalSubspace;
    globalToGridLocal = std::make_unique<Acts::GridAccess::Affine3Transformed<
        Acts::GridAccess::GlobalSubspace<Args...>>>(std::move(globalSubspace),
                                                    transform);
  } else {
    globalToGridLocal =
        std::make_unique<Acts::GridAccess::GlobalSubspace<Args...>>();
  }
  return globalToGridLocal;
}

}  // namespace

nlohmann::json Acts::GridAccessJsonConverter::toJson(
    const GridAccess::IGlobalToGridLocal& globalToGridLocal) {
  nlohmann::json jGlobalToGridLocal;

  std::array<bool, 2u> transformOptions = {false, true};

  // One dimensional sub spaces
  const std::tuple<
      GridAccess::GlobalSubspace<binX>, GridAccess::GlobalSubspace<binY>,
      GridAccess::GlobalSubspace<binZ>, GridAccess::GlobalSubspace<binR>,
      GridAccess::GlobalSubspace<binPhi>, GridAccess::GlobalSubspace<binEta>>
      oneDimSubspaces = {};

  for (bool transform : transformOptions) {
    encodeSubspaces(jGlobalToGridLocal, globalToGridLocal, transform,
                    oneDimSubspaces);
    if (!jGlobalToGridLocal.empty()) {
      return jGlobalToGridLocal;
    }
  }

  // Useful two dimensional sub spaces
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
    encodeSubspaces(jGlobalToGridLocal, globalToGridLocal, transform,
                    twoDimSubspaces);
    if (!jGlobalToGridLocal.empty()) {
      return jGlobalToGridLocal;
    }
  }
  return jGlobalToGridLocal;
}

std::unique_ptr<Acts::GridAccess::IGlobalToGridLocal>
Acts::GridAccessJsonConverter::globalToGridLocalFromJson(
    const nlohmann::json& jGlobalToGridLocal) {
  std::unique_ptr<Acts::GridAccess::IGlobalToGridLocal> globalToGridLocal =
      nullptr;

  std::vector<BinningValue> accessors =
      jGlobalToGridLocal.at("accessors").get<std::vector<BinningValue>>();

  // Switch and fill for 1D
  if (accessors.size() == 1u) {
    switch (accessors[0]) {
      case binX:
        globalToGridLocal = decodeSubspace<binX>(jGlobalToGridLocal);
        break;
      case binY:
        globalToGridLocal = decodeSubspace<binY>(jGlobalToGridLocal);
        break;
      case binZ:
        globalToGridLocal = decodeSubspace<binZ>(jGlobalToGridLocal);
        break;
      case binR:
        globalToGridLocal = decodeSubspace<binR>(jGlobalToGridLocal);
        break;
      case binPhi:
        globalToGridLocal = decodeSubspace<binPhi>(jGlobalToGridLocal);
        break;
      case binEta:
        globalToGridLocal = decodeSubspace<binEta>(jGlobalToGridLocal);
        break;
      default:
        // globalToGridLocal = nullptr;
        break;
    }
  }

  // Switch and fill for 2D
  if (accessors.size() == 2u) {
    if (accessors == std::vector<BinningValue>{binX, binY}) {
      globalToGridLocal = decodeSubspace<binX, binY>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binY, binX}) {
      globalToGridLocal = decodeSubspace<binY, binX>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binX, binZ}) {
      globalToGridLocal = decodeSubspace<binX, binZ>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binZ, binX}) {
      globalToGridLocal = decodeSubspace<binZ, binX>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binY, binZ}) {
      globalToGridLocal = decodeSubspace<binY, binZ>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binZ, binY}) {
      globalToGridLocal = decodeSubspace<binZ, binY>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binR, binPhi}) {
      globalToGridLocal = decodeSubspace<binR, binPhi>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binPhi, binR}) {
      globalToGridLocal = decodeSubspace<binPhi, binR>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binZ, binPhi}) {
      globalToGridLocal = decodeSubspace<binZ, binPhi>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binPhi, binZ}) {
      globalToGridLocal = decodeSubspace<binPhi, binZ>(jGlobalToGridLocal);
    }
    // else globalToGridLocal = nullptr;
  }
  return globalToGridLocal;
}

nlohmann::json Acts::GridAccessJsonConverter::toJson(
    const GridAccess::IBoundToGridLocal& boundToGridLocal) {
  nlohmann::json jBoundtoGridLocal;

  auto localSubSpace0 =
      dynamic_cast<const GridAccess::LocalSubspace<0u>*>(&boundToGridLocal);
  if (localSubSpace0 != nullptr) {
    jBoundtoGridLocal["type"] = "subspace";
    jBoundtoGridLocal["accessors"] = localSubSpace0->accessors;
  }

  auto localSubSpace1 =
      dynamic_cast<const GridAccess::LocalSubspace<1u>*>(&boundToGridLocal);
  if (localSubSpace1 != nullptr) {
    jBoundtoGridLocal["type"] = "subspace";
    jBoundtoGridLocal["accessors"] = localSubSpace1->accessors;
  }

  auto localSubSpace01 =
      dynamic_cast<const GridAccess::LocalSubspace<0u, 1u>*>(&boundToGridLocal);
  if (localSubSpace01 != nullptr) {
    jBoundtoGridLocal["type"] = "subspace";
    jBoundtoGridLocal["accessors"] = localSubSpace01->accessors;
  }

  auto localSubSpace10 =
      dynamic_cast<const GridAccess::LocalSubspace<1u, 0u>*>(&boundToGridLocal);
  if (localSubSpace10 != nullptr) {
    jBoundtoGridLocal["type"] = "subspace";
    jBoundtoGridLocal["accessors"] = localSubSpace10->accessors;
  }

  auto boundCylinderToZPhi =
      dynamic_cast<const GridAccess::BoundCylinderToZPhi*>(&boundToGridLocal);
  if (boundCylinderToZPhi != nullptr) {
    jBoundtoGridLocal["type"] = "cylinder_to_zphi";
    jBoundtoGridLocal["radius"] = boundCylinderToZPhi->radius;
    jBoundtoGridLocal["shift"] = boundCylinderToZPhi->shift;
  }

  return jBoundtoGridLocal;
}

std::unique_ptr<Acts::GridAccess::IBoundToGridLocal>
Acts::GridAccessJsonConverter::boundToGridLocalFromJson(
    const nlohmann::json& jBoundtoGridLocal) {
  std::unique_ptr<Acts::GridAccess::IBoundToGridLocal> boundToGridLocal =
      nullptr;
  std::string type = jBoundtoGridLocal.at("type").get<std::string>();
  if (type == "subspace") {
    std::vector<std::size_t> accessors =
        jBoundtoGridLocal.at("accessors").get<std::vector<std::size_t>>();
    if (accessors.size() == 1 && accessors[0] == 0) {
      boundToGridLocal =
          std::make_unique<Acts::GridAccess::LocalSubspace<0u>>();
    } else if (accessors.size() == 1 && accessors[0] == 1) {
      boundToGridLocal =
          std::make_unique<Acts::GridAccess::LocalSubspace<1u>>();
    } else if (accessors.size() == 2 && accessors[0] == 0 &&
               accessors[1] == 1) {
      boundToGridLocal =
          std::make_unique<Acts::GridAccess::LocalSubspace<0u, 1u>>();
    } else if (accessors.size() == 2 && accessors[0] == 1 &&
               accessors[1] == 0) {
      boundToGridLocal =
          std::make_unique<Acts::GridAccess::LocalSubspace<1u, 0u>>();
    }
  } else if (type == "cylinder_to_zphi") {
    ActsScalar radius = jBoundtoGridLocal.at("radius").get<ActsScalar>();
    ActsScalar shift = jBoundtoGridLocal.at("shift").get<ActsScalar>();
    boundToGridLocal =
        std::make_unique<Acts::GridAccess::BoundCylinderToZPhi>(radius, shift);
  }
  return boundToGridLocal;
}

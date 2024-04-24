// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/GridJsonConverter.hpp"

#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Plugins/Json/UtilitiesJsonConverter.hpp"
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
    const Subspace& /*subspace*/) {
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
std::unique_ptr<const Acts::GridAccess::GlobalSubspace<Args...>> decodeSubspace(
    const nlohmann::json& /*j*/) {
  return std::make_unique<const Acts::GridAccess::GlobalSubspace<Args...>>();
}

template <Acts::BinningValue... Args>
std::unique_ptr<const Acts::GridAccess::Affine3Transformed<
    Acts::GridAccess::GlobalSubspace<Args...>>>
decodeTransformedSubspace(const nlohmann::json& jGlobalToGridLocal) {
  Acts::Transform3 transform = Acts::Transform3JsonConverter::fromJson(
      jGlobalToGridLocal.at("transform"));
  Acts::GridAccess::GlobalSubspace<Args...> globalSubspace;
  return std::make_unique<const Acts::GridAccess::Affine3Transformed<
      Acts::GridAccess::GlobalSubspace<Args...>>>(std::move(globalSubspace),
                                                  transform);
}

template <Acts::BinningValue... Args>
std::unique_ptr<const Acts::GridAccess::IGlobalToGridLocal>
decodeGeneralSubspace(const nlohmann::json& jGlobalToGridLocal) {
  if (jGlobalToGridLocal.find("transform") != jGlobalToGridLocal.end()) {
    return decodeTransformedSubspace<Args...>(jGlobalToGridLocal);
  }
  return decodeSubspace<Args...>(jGlobalToGridLocal);
}

template <typename Delegate, Acts::BinningValue... Args>
void decorateGlobalDelegate(Delegate& delegate,
                            const nlohmann::json& jGlobalToGridLocal) {
  // The delegate has already been connected
  if (delegate.connected()) {
    return;
  }

  // Get the transform for json
  bool hasTransform =
      jGlobalToGridLocal.find("transform") != jGlobalToGridLocal.end();

  // Get the accessors
  std::vector<Acts::BinningValue> accessors =
      jGlobalToGridLocal.at("accessors").get<std::vector<Acts::BinningValue>>();

  // One dimensional setting
  if constexpr (sizeof...(Args) == 1u) {
    if (std::get<0>(std::forward_as_tuple(Args...)) == accessors[0]) {
      if (hasTransform) {
        using TransformedSubspace = Acts::GridAccess::Affine3Transformed<
            Acts::GridAccess::GlobalSubspace<Args...>>;
        auto globalToGridLocal =
            decodeTransformedSubspace<Args...>(jGlobalToGridLocal);
        delegate.template connect<&TransformedSubspace::toGridLocal>(
            std::move(globalToGridLocal));
      } else {
        auto globalToGridLocal = decodeSubspace<Args...>(jGlobalToGridLocal);
        delegate.template connect<
            &Acts::GridAccess::GlobalSubspace<Args...>::toGridLocal>(
            std::move(globalToGridLocal));
      }
    }
  }

  // Two-dimensional setting
  if constexpr (sizeof...(Args) == 2u) {
    if (std::get<0>(std::forward_as_tuple(Args...)) == accessors[0] &&
        std::get<1>(std::forward_as_tuple(Args...)) == accessors[1]) {
      if (hasTransform) {
        using TransformedSubspace = Acts::GridAccess::Affine3Transformed<
            Acts::GridAccess::GlobalSubspace<Args...>>;
        auto globalToGridLocal =
            decodeTransformedSubspace<Args...>(jGlobalToGridLocal);
        delegate.template connect<&TransformedSubspace::toGridLocal>(
            std::move(globalToGridLocal));
      } else {
        auto globalToGridLocal = decodeSubspace<Args...>(jGlobalToGridLocal);
        delegate.template connect<
            &Acts::GridAccess::GlobalSubspace<Args...>::toGridLocal>(
            std::move(globalToGridLocal));
      }
    }
  }
}

template <Acts::BinningValue... Args>
void decorateGlobal1DimDelegate(
    Acts::GridAccess::GlobalToGridLocal1DimDelegate& delegate,
    const nlohmann::json& jGlobalToGridLocal) {
  (((decorateGlobalDelegate<Acts::GridAccess::GlobalToGridLocal1DimDelegate,
                            Args>(delegate, jGlobalToGridLocal))),
   ...);
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

std::unique_ptr<const Acts::GridAccess::IGlobalToGridLocal>
Acts::GridAccessJsonConverter::globalToGridLocalFromJson(
    const nlohmann::json& jGlobalToGridLocal) {
  std::unique_ptr<const Acts::GridAccess::IGlobalToGridLocal>
      globalToGridLocal = nullptr;

  std::vector<BinningValue> accessors =
      jGlobalToGridLocal.at("accessors").get<std::vector<BinningValue>>();

  // Switch and fill for 1D
  if (accessors.size() == 1u) {
    switch (accessors[0]) {
      case binX:
        globalToGridLocal = decodeGeneralSubspace<binX>(jGlobalToGridLocal);
        break;
      case binY:
        globalToGridLocal = decodeGeneralSubspace<binY>(jGlobalToGridLocal);
        break;
      case binZ:
        globalToGridLocal = decodeGeneralSubspace<binZ>(jGlobalToGridLocal);
        break;
      case binR:
        globalToGridLocal = decodeGeneralSubspace<binR>(jGlobalToGridLocal);
        break;
      case binPhi:
        globalToGridLocal = decodeGeneralSubspace<binPhi>(jGlobalToGridLocal);
        break;
      case binEta:
        globalToGridLocal = decodeGeneralSubspace<binEta>(jGlobalToGridLocal);
        break;
      default:
        // globalToGridLocal = nullptr;
        break;
    }
  }

  // Switch and fill for 2D
  if (accessors.size() == 2u) {
    if (accessors == std::vector<BinningValue>{binX, binY}) {
      globalToGridLocal = decodeGeneralSubspace<binX, binY>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binY, binX}) {
      globalToGridLocal = decodeGeneralSubspace<binY, binX>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binX, binZ}) {
      globalToGridLocal = decodeGeneralSubspace<binX, binZ>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binZ, binX}) {
      globalToGridLocal = decodeGeneralSubspace<binZ, binX>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binY, binZ}) {
      globalToGridLocal = decodeGeneralSubspace<binY, binZ>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binZ, binY}) {
      globalToGridLocal = decodeGeneralSubspace<binZ, binY>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binR, binPhi}) {
      globalToGridLocal =
          decodeGeneralSubspace<binR, binPhi>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binPhi, binR}) {
      globalToGridLocal =
          decodeGeneralSubspace<binPhi, binR>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binZ, binPhi}) {
      globalToGridLocal =
          decodeGeneralSubspace<binZ, binPhi>(jGlobalToGridLocal);
    } else if (accessors == std::vector<BinningValue>{binPhi, binZ}) {
      globalToGridLocal =
          decodeGeneralSubspace<binPhi, binZ>(jGlobalToGridLocal);
    }
    // else globalToGridLocal = nullptr;
  }
  return globalToGridLocal;
}

Acts::GridAccess::GlobalToGridLocal1DimDelegate
Acts::GridAccessJsonConverter::globalToGridLocal1DimDelegateFromJson(
    const nlohmann::json& jGlobalToGridLocal) {
  // Peek into json to check the right dimension
  if (jGlobalToGridLocal.at("accessors").size() != 1u) {
    throw std::invalid_argument(
        "GridAccessJsonConverter: json input does not describe 1D case.");
  }
  // Unroll the decoration
  Acts::GridAccess::GlobalToGridLocal1DimDelegate delegate;
  decorateGlobal1DimDelegate<binX, binY, binZ, binR, binPhi, binEta>(
      delegate, jGlobalToGridLocal);
  return delegate;
}

Acts::GridAccess::GlobalToGridLocal2DimDelegate
Acts::GridAccessJsonConverter::globalToGridLocal2DimDelegateFromJson(
    const nlohmann::json& jGlobalToGridLocal) {
  // Peek into json to check the right dimension
  if (jGlobalToGridLocal.at("accessors").size() != 2u) {
    throw std::invalid_argument(
        "GridAccessJsonConverter: json input does not describe 2D case.");
  }
  // Unroll the decoration
  Acts::GridAccess::GlobalToGridLocal2DimDelegate delegate;
  // Only the matching one will be applied, matching condition is checked inside
  // the call - may unroll this es well
  decorateGlobalDelegate<Acts::GridAccess::GlobalToGridLocal2DimDelegate, binX,
                         binY>(delegate, jGlobalToGridLocal);
  decorateGlobalDelegate<Acts::GridAccess::GlobalToGridLocal2DimDelegate, binY,
                         binX>(delegate, jGlobalToGridLocal);
  decorateGlobalDelegate<Acts::GridAccess::GlobalToGridLocal2DimDelegate, binX,
                         binZ>(delegate, jGlobalToGridLocal);
  decorateGlobalDelegate<Acts::GridAccess::GlobalToGridLocal2DimDelegate, binZ,
                         binX>(delegate, jGlobalToGridLocal);
  decorateGlobalDelegate<Acts::GridAccess::GlobalToGridLocal2DimDelegate, binY,
                         binZ>(delegate, jGlobalToGridLocal);
  decorateGlobalDelegate<Acts::GridAccess::GlobalToGridLocal2DimDelegate, binZ,
                         binY>(delegate, jGlobalToGridLocal);
  decorateGlobalDelegate<Acts::GridAccess::GlobalToGridLocal2DimDelegate, binR,
                         binPhi>(delegate, jGlobalToGridLocal);
  decorateGlobalDelegate<Acts::GridAccess::GlobalToGridLocal2DimDelegate,
                         binPhi, binR>(delegate, jGlobalToGridLocal);
  decorateGlobalDelegate<Acts::GridAccess::GlobalToGridLocal2DimDelegate, binZ,
                         binPhi>(delegate, jGlobalToGridLocal);
  decorateGlobalDelegate<Acts::GridAccess::GlobalToGridLocal2DimDelegate,
                         binPhi, binZ>(delegate, jGlobalToGridLocal);
  return delegate;
}

nlohmann::json Acts::GridAccessJsonConverter::toJson(
    const GridAccess::IBoundToGridLocal& boundToGridLocal) {
  nlohmann::json jBoundToGridLocal;

  auto localSubSpace0 =
      dynamic_cast<const GridAccess::LocalSubspace<0u>*>(&boundToGridLocal);
  if (localSubSpace0 != nullptr) {
    jBoundToGridLocal["type"] = "subspace";
    jBoundToGridLocal["accessors"] = localSubSpace0->accessors;
  }

  auto localSubSpace1 =
      dynamic_cast<const GridAccess::LocalSubspace<1u>*>(&boundToGridLocal);
  if (localSubSpace1 != nullptr) {
    jBoundToGridLocal["type"] = "subspace";
    jBoundToGridLocal["accessors"] = localSubSpace1->accessors;
  }

  auto localSubSpace01 =
      dynamic_cast<const GridAccess::LocalSubspace<0u, 1u>*>(&boundToGridLocal);
  if (localSubSpace01 != nullptr) {
    jBoundToGridLocal["type"] = "subspace";
    jBoundToGridLocal["accessors"] = localSubSpace01->accessors;
  }

  auto localSubSpace10 =
      dynamic_cast<const GridAccess::LocalSubspace<1u, 0u>*>(&boundToGridLocal);
  if (localSubSpace10 != nullptr) {
    jBoundToGridLocal["type"] = "subspace";
    jBoundToGridLocal["accessors"] = localSubSpace10->accessors;
  }

  auto boundCylinderToZPhi =
      dynamic_cast<const GridAccess::BoundCylinderToZPhi*>(&boundToGridLocal);
  if (boundCylinderToZPhi != nullptr) {
    jBoundToGridLocal["type"] = "cylinder_to_zphi";
    jBoundToGridLocal["radius"] = boundCylinderToZPhi->radius;
    jBoundToGridLocal["shift"] = boundCylinderToZPhi->shift;
  }

  return jBoundToGridLocal;
}

std::unique_ptr<Acts::GridAccess::IBoundToGridLocal>
Acts::GridAccessJsonConverter::boundToGridLocalFromJson(
    const nlohmann::json& jBoundToGridLocal) {
  std::unique_ptr<Acts::GridAccess::IBoundToGridLocal> boundToGridLocal =
      nullptr;
  std::string type = jBoundToGridLocal.at("type").get<std::string>();
  if (type == "subspace") {
    std::vector<std::size_t> accessors =
        jBoundToGridLocal.at("accessors").get<std::vector<std::size_t>>();
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
    ActsScalar radius = jBoundToGridLocal.at("radius").get<ActsScalar>();
    ActsScalar shift = jBoundToGridLocal.at("shift").get<ActsScalar>();
    boundToGridLocal =
        std::make_unique<Acts::GridAccess::BoundCylinderToZPhi>(radius, shift);
  }
  return boundToGridLocal;
}

Acts::GridAccess::BoundToGridLocal1DimDelegate
Acts::GridAccessJsonConverter::boundToGridLocal1DimDelegateFromJson(
    const nlohmann::json& jBoundToGridLocal) {
  Acts::GridAccess::BoundToGridLocal1DimDelegate delegate;

  std::string type = jBoundToGridLocal.at("type").get<std::string>();
  if (type == "subspace") {
    std::vector<std::size_t> accessors =
        jBoundToGridLocal.at("accessors").get<std::vector<std::size_t>>();
    // Safety check
    if (accessors.size() != 1u) {
      throw std::invalid_argument(
          "GridAccessJsonConverter: json input does not describe 1D case.");
    }
    // Specify the type
    if (accessors[0] == 0) {
      auto boundToGridLocal =
          std::make_unique<const Acts::GridAccess::LocalSubspace<0u>>();
      delegate.connect<&Acts::GridAccess::LocalSubspace<0u>::toGridLocal>(
          std::move(boundToGridLocal));
    } else if (accessors[0] == 1) {
      auto boundToGridLocal =
          std::make_unique<const Acts::GridAccess::LocalSubspace<1u>>();
      delegate.connect<&Acts::GridAccess::LocalSubspace<1u>::toGridLocal>(
          std::move(boundToGridLocal));
    }
  }
  return delegate;
}

Acts::GridAccess::BoundToGridLocal2DimDelegate
Acts::GridAccessJsonConverter::boundToGridLocal2DimDelegateFromJson(
    const nlohmann::json& jBoundToGridLocal) {
  Acts::GridAccess::BoundToGridLocal2DimDelegate delegate;

  std::string type = jBoundToGridLocal.at("type").get<std::string>();
  if (type == "subspace") {
    std::vector<std::size_t> accessors =
        jBoundToGridLocal.at("accessors").get<std::vector<std::size_t>>();

    // Safety check
    if (accessors.size() != 2u) {
      throw std::invalid_argument(
          "GridAccessJsonConverter: json input does not describe 2D case.");
    }
    if (accessors[0] == 0u && accessors[1] == 1u) {
      auto boundToGridLocal =
          std::make_unique<const Acts::GridAccess::LocalSubspace<0u, 1u>>();
      delegate.connect<&Acts::GridAccess::LocalSubspace<0u, 1u>::toGridLocal>(
          std::move(boundToGridLocal));
    } else if (accessors[0] == 1u && accessors[1] == 0u) {
      auto boundToGridLocal =
          std::make_unique<const Acts::GridAccess::LocalSubspace<1u, 0u>>();
      delegate.connect<&Acts::GridAccess::LocalSubspace<1u, 0u>::toGridLocal>(
          std::move(boundToGridLocal));
    }
  } else if (type == "cylinder_to_zphi") {
    ActsScalar radius = jBoundToGridLocal.at("radius").get<ActsScalar>();
    ActsScalar shift = jBoundToGridLocal.at("shift").get<ActsScalar>();
    auto boundToGridLocal =
        std::make_unique<const Acts::GridAccess::BoundCylinderToZPhi>(radius,
                                                                      shift);
    delegate.connect<&Acts::GridAccess::BoundCylinderToZPhi::toGridLocal>(
        std::move(boundToGridLocal));
  }

  return delegate;
}

// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/GridJsonConverter.hpp"

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

template <typename Subspace>
void encodeSubspace(nlohmann::json& j,
                    const Acts::GridAccess::IGlobalToGridLocal& ga,
                    const Subspace& /*unused*/) {
  const Subspace* subspace = dynamic_cast<const Subspace*>(&ga);
  if (subspace != nullptr) {
    j["type"] = "subspace";
    j["accessors"] = subspace->accessors;
  }
}

template <typename... Args>
void encodeSubspaces(nlohmann::json& j,
                     const Acts::GridAccess::IGlobalToGridLocal& ga,
                     std::tuple<Args...>& t) {
  std::apply([&](auto&&... vals) { (encodeSubspace(j, ga, vals), ...); }, t);                    
}

nlohmann::json Acts::GridAccessJsonConverter::toJson(
    const GridAccess::IGlobalToGridLocal& ga) {
  nlohmann::json jga;

  // One dimensional sub spaces
  const std::tuple<
      GridAccess::GlobalSubspace<binX>, GridAccess::GlobalSubspace<binY>,
      GridAccess::GlobalSubspace<binZ>, GridAccess::GlobalSubspace<binR>,
      GridAccess::GlobalSubspace<binPhi>, GridAccess::GlobalSubspace<binEta>>
      oneDimAccessors = {};

  encodeSubspaces(jga, ga, oneDimAccessors);
  if (!jga.empty()) {
    return jga;
  }

  // Usfeul two dimensional sub spaces 
  const std::tuple<
    GridAccess::GlobalSubspace<binX, binY>, GridAccess::GlobalSubspace<binY, binX>,
    GridAccess::GlobalSubspace<binX, binZ>, GridAccess::GlobalSubspace<binZ, binX>,
    GridAccess::GlobalSubspace<binY, binZ>, GridAccess::GlobalSubspace<binZ, binY>,
    GridAccess::GlobalSubspace<binR, binPhi>, GridAccess::GlobalSubspace<binPhi, binR>, 
    GridAccess::GlobalSubspace<binZ, binPhi>, GridAccess::GlobalSubspace<binPhi, binZ>>
    twoDimAccessors = {};

  encodeSubspaces(jga, ga, twoDimAccessors);
  if (!jga.empty()) {
    return jga;
  }



  return jga;
}

std::unique_ptr<Acts::GridAccess::IGlobalToGridLocal>
Acts::GridAccessJsonConverter::globalToGridLocalFromJson(
    const nlohmann::json& jga) {
  std::unique_ptr<Acts::GridAccess::IGlobalToGridLocal> ga;
  return ga;
}

nlohmann::json Acts::GridAccessJsonConverter::toJson(
    const GridAccess::IBoundToGridLocal& la) {
  nlohmann::json jla;

  auto localSubSpace0 = dynamic_cast<const GridAccess::LocalSubspace<0u>*>(&la);
  if (localSubSpace0 != nullptr) {
    jla["type"] = "subspace";
    jla["accessors"] = localSubSpace0->accessors;
  }

  auto localSubSpace1 = dynamic_cast<const GridAccess::LocalSubspace<1u>*>(&la);
  if (localSubSpace1 != nullptr) {
    jla["type"] = "subspace";
    jla["accessors"] = localSubSpace1->accessors;
  }

  auto localSubSpace01 =
      dynamic_cast<const GridAccess::LocalSubspace<0u, 1u>*>(&la);
  if (localSubSpace01 != nullptr) {
    jla["type"] = "subspace";
    jla["accessors"] = localSubSpace01->accessors;
  }

  auto localSubSpace10 =
      dynamic_cast<const GridAccess::LocalSubspace<1u, 0u>*>(&la);
  if (localSubSpace10 != nullptr) {
    jla["type"] = "subspace";
    jla["accessors"] = localSubSpace10->accessors;
  }

  auto boundCylinderToZPhi =
      dynamic_cast<const GridAccess::BoundCylinderToZPhi*>(&la);
  if (boundCylinderToZPhi != nullptr) {
    jla["type"] = "cylinder_to_zphi";
    jla["radius"] = boundCylinderToZPhi->radius;
    jla["shift"] = boundCylinderToZPhi->shift;
  }

  return jla;
}

std::unique_ptr<Acts::GridAccess::IBoundToGridLocal>
Acts::GridAccessJsonConverter::boundToGridLocalFromJson(
    const nlohmann::json& jla) {
  std::unique_ptr<Acts::GridAccess::IBoundToGridLocal> la;
  std::string type = jla.at("type").get<std::string>();
  if (type == "subspace") {
    std::vector<std::size_t> accessors =
        jla.at("accessors").get<std::vector<std::size_t>>();
    if (accessors.size() == 1 && accessors[0] == 0) {
      la = std::make_unique<Acts::GridAccess::LocalSubspace<0u>>();
    } else if (accessors.size() == 1 && accessors[0] == 1) {
      la = std::make_unique<Acts::GridAccess::LocalSubspace<1u>>();
    } else if (accessors.size() == 2 && accessors[0] == 0 &&
               accessors[1] == 1) {
      la = std::make_unique<Acts::GridAccess::LocalSubspace<0u, 1u>>();
    } else if (accessors.size() == 2 && accessors[0] == 1 &&
               accessors[1] == 0) {
      la = std::make_unique<Acts::GridAccess::LocalSubspace<1u, 0u>>();
    }
  } else if (type == "cylinder_to_zphi") {
    ActsScalar radius = jla.at("radius").get<ActsScalar>();
    ActsScalar shift = jla.at("shift").get<ActsScalar>();
    la = std::make_unique<Acts::GridAccess::BoundCylinderToZPhi>(radius, shift);
  }
  return la;
}

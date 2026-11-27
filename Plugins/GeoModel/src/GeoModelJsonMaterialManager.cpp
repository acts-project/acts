// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/GeoModel/GeoModelJsonMaterialManager.hpp"

#include "ActsPlugins/Json/ActsJson.hpp"

#include <fstream>
#include <iostream>

namespace ActsPlugins {
GeoModelJsonMaterialManager* GeoModelJsonMaterialManager::getManager() {
  if (s_instance == nullptr) {
    s_instance = new GeoModelJsonMaterialManager();
  }
  return dynamic_cast<GeoModelJsonMaterialManager*>(s_instance);
}

bool GeoModelJsonMaterialManager::loadMaterialMap(const std::string& filePath) {
  std::ifstream istr{filePath};
  if (!istr.good()) {
    ACTS_ERROR("Failed to open '" << filePath << "'");
    return false;
  }
  nlohmann::json payload{};
  istr >> payload;

  /// Load first the periodic table
  for (const auto& elementPayload : payload["periodicTable"]) {
    const std::string eleName = elementPayload["name"];
    const std::string symbol = elementPayload["symbol"];
    const unsigned Z = elementPayload["charge"];
    const double A = elementPayload["atomicNumber"];

    if (isElementDefined(eleName)) {
      ACTS_ERROR("Cannot redefine the element '" << eleName << "'");
      return false;
    }
    ACTS_DEBUG("Add new chemical element " << eleName << "(" << symbol
                                           << ") with Z=" << Z
                                           << " and A=" << A);
    addElement(eleName, symbol, Z, A);
  }
  ACTS_DEBUG(
      "\n\n\n Elements are all defined now. Start parsing the materials");
  setMaterialNamespace("std");
  /// Then load the matieral table
  for (const auto& materialPayload : payload["materials"]) {
    const std::string name = materialPayload["name"];
    const double density = materialPayload["rho"];
    if (materialPayload.find("namespace") != materialPayload.end()) {
      const std::string nSpace = materialPayload["namespace"];
      setMaterialNamespace(nSpace);
    }
    addMaterial(name, density);
    ACTS_DEBUG("Define new material " << name << " with density " << density
                                      << " in namespace "
                                      << materialNameSpace());
    for (const auto& compPayload : materialPayload["composition"]) {
      const std::string component = compPayload["component"];
      const double parts = compPayload["parts"];
      ACTS_VERBOSE("    Add " << component << " as material component with "
                              << parts);
      addMatComponent(component, parts);
    }
    lockMaterial();
    setMaterialNamespace("std");
  }
  return true;
}
}  // namespace ActsPlugins

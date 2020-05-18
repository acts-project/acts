// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Geometry/HierarchicalGeometryContainer.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <nlohmann/json.hpp>

#include <map>

namespace Acts {

/// Convert Hierarchical Object Container to Json file and vice-versa
template <typename object_t>
struct JsonHierarchicalGeometryContainerConverter {
 public:
  /// The object key
  std::string datakey = "Object";

  /// The loglevel
  Acts::Logging::Level logLevel = Acts::Logging::Level::INFO;

  /// Convert json map to Hierarchical Object Container
  ///
  /// @param map The indexed Hierarchical Object in json format
  /// @param fromJson Function that return Hierarchical Object from json
  HierarchicalGeometryContainer<object_t> jsonToHierarchicalContainer(
      const nlohmann::json& map,
      std::function<object_t(const Acts::GeometryID&, const nlohmann::json&)>
          fromJson) const;

  /// Convert json map to Hierarchical Object Container
  ///
  /// @param hObject The Hierarchical Object Container
  /// @param toJson Function that convert Hierarchical Object to json
  nlohmann::json hierarchicalObjectToJson(
      const HierarchicalGeometryContainer<object_t>& hObject,
      std::function<nlohmann::json(const object_t&)> toJson) const;

 private:
  /// The detector tag
  static constexpr char m_detkey[] = "detector";
  /// The volume identification string
  static constexpr char m_volkey[] = "volumes";
  /// The boundary surface string
  static constexpr char m_boukey[] = "boundaries";
  /// The layer identification string
  static constexpr char m_laykey[] = "layers";
  /// The approach identification string
  static constexpr char m_appkey[] = "approach";
  /// The sensitive identification string
  static constexpr char m_senkey[] = "sensitive";
  /// The representing idntification string
  static constexpr char m_repkey[] = "representing";
  /// The name identification
  static constexpr char m_namekey[] = "Name";

  // helper function to create geometry id from json keys
  Acts::GeometryID makeId(std::string volume, std::string boundary,
                          std::string layer, std::string approach,
                          std::string sensitive) const {
    return GeometryID()
        .setVolume(std::stoi(volume))
        .setBoundary(std::stoi(boundary))
        .setLayer(std::stoi(layer))
        .setApproach(std::stoi(approach))
        .setSensitive(std::stoi(sensitive));
  }
};

}  // namespace Acts
#include "Acts/Plugins/Json/JsonHierarchicalGeometryContainerConverter.ipp"

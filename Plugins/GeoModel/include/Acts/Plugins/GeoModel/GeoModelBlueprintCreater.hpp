// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/Blueprint.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace Acts {

struct GeoModelTree;

class GeoModelBlueprintCreater {
 public:
  /// The configuration struct
  struct Config {};

  /// The Options struct
  struct Options {
    /// The table name of the blueprint in the database
    std::string table;
    /// The top level node name
    std::string topEntry;
    /// Optionally override the top node bounds
    std::vector<ActsScalar> topNodeBounds;
  };

  /// The Blueprint return object
  struct Blueprint {
    std::string name;
    std::unique_ptr<Acts::Experimental::Blueprint::Node> topNode;

    /// Access to the top node
    const Acts::Experimental::Blueprint::Node& node() const {
      if (topNode == nullptr) {
        throw std::runtime_error(
            "GeoModelBlueprintCreater::Blueprint: No top node created");
      }
      return *(topNode.get());
    }
  };

  // Blueprint table is:
  // cursor.execute("CREATE TABLE Blueprint(id INT, type TEXT, name TEXT, bounds
  // TEXT, internals TEXT, binnings TEXT, materials TEXT)")
  //
  struct TableEntry {
    int id;
    std::string type;
    std::string name;
    std::vector<std::string> bounds;
    std::vector<std::string> internals;
    std::vector<std::string> binnings;
    std::vector<std::string> materials;
  };

  /// The GeoModel blueprint creater from the database
  ///
  /// @param cfg the configuration struct
  /// @param mlogger a screen output logger
  GeoModelBlueprintCreater(const Config& cfg,
                           std::unique_ptr<const Logger> mlogger =
                               getDefaultLogger("GeoModelBlueprintCreater",
                                                Acts::Logging::INFO));

  /// Method that reads the GeoModel blueprint from database
  ///
  /// @param gctx the geometry context
  /// @param gmTree the GeoModel tree
  /// @param options the options
  Blueprint create(const GeometryContext& gctx, const GeoModelTree& gmTree,
                   const Options& options) const;

 private:
  /// Create a blueprint node from a table entry
  ///
  /// @param entry the table entry
  /// @param tableEntryMap the map of table entries allows construction of children
  /// @param motherExtent the raw extent of the mother
  ///
  /// @return a newly created node
  std::unique_ptr<Experimental::Blueprint::Node> createNode(
      const TableEntry& entry,
      const std::map<std::string, TableEntry>& tableEntryMap,
      const Extent& motherExtent = Extent()) const;

  /// @brief Parse bound value string from the database
  ///
  /// @param boundsEntry in the database
  /// @param rawValuesMother the raw bound values of the mother
  /// @param motherExtent the raw extent of the mother
  ///
  /// @return The bounds type, raw bound values, deduced bound values, and a translation vector
  std::tuple<VolumeBounds::BoundsType, Extent, std::vector<ActsScalar>, Vector3>
  parseBounds(const std::vector<std::string>& boundsEntry,
              const Extent& motherExtent) const;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to the logger
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

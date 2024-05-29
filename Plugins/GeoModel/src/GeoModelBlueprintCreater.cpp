// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelBlueprintCreater.hpp"

#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Enumerate.hpp"

namespace {

/// @brief Helper to transform binning string to BinningValue enum
Acts::BinningValue toBinningValue(const std::string& binning) {
  if (binning == "x") {
    return Acts::BinningValue::binX;
  } else if (binning == "y") {
    return Acts::BinningValue::binY;
  } else if (binning == "z") {
    return Acts::BinningValue::binZ;
  } else if (binning == "r") {
    return Acts::BinningValue::binR;
  } else if (binning == "phi") {
    return Acts::BinningValue::binPhi;
  } else {
    throw std::invalid_argument(
        "GeoModelBlueprintCreater: Unknown binning value '" + binning + "'");
  }
}

/// @brief Split a string into a vector of strings
///
/// @param entry the string to be split
/// @param deliminater split indicator
///
/// @return a vector of strings
std::vector<std::string> splitString(const std::string& entry,
                                     const std::string& deliminater) {
  std::vector<std::string> result;
  std::string currentEntry = entry;
  size_t pos = 0;
  while ((pos = currentEntry.find(deliminater)) != std::string::npos) {
    auto found = currentEntry.substr(0, pos);
    if (!found.empty()) {
      result.push_back(currentEntry.substr(0, pos));
    }
    currentEntry.erase(0, pos + deliminater.length());
  }
  result.push_back(currentEntry);
  return result;
}

}  // namespace

Acts::GeoModelBlueprintCreater::GeoModelBlueprintCreater(
    const Config&, std::unique_ptr<const Logger> mlogger)
    : m_logger(std::move(mlogger)) {}

Acts::GeoModelBlueprintCreater::Blueprint
Acts::GeoModelBlueprintCreater::create(const GeometryContext& gctx,
                                       const GeoModelTree& gmTree,
                                       const Options& options) const {
  // The blueprint to be created
  Acts::GeoModelBlueprintCreater::Blueprint blueprint;

  // The GeoModel tree must have a reader
  if (gmTree.geoReader == nullptr) {
    throw std::invalid_argument(
        "GeoModelBlueprintCreater: GeoModelTree has no GeoModelReader");
  }

  auto blueprintTable = gmTree.geoReader->getTableFromTableName(options.table);

  std::map<std::string, TableEntry> blueprintTableMap;

  // In order to guarantee the correct order
  // of the building sequence we need to build a representative tree first
  for (const auto& line : blueprintTable) {
    if (line.size() != 7u) {
      throw std::invalid_argument(
          "GeoModelBlueprintCreater: Blueprint table has wrong number of "
          "columns");
    }

    int volumeId = std::stoi(line.at(0));
    std::string volumeType = line.at(1);
    std::string volumeName = line.at(2);
    // The bit more complicated strings from the database
    std::vector<std::string> volumeBounds = splitString(line.at(3), ";");
    std::vector<std::string> volumeInternals = splitString(line.at(4), ";");
    std::vector<std::string> volumeBinnings = splitString(line.at(5), ";");
    std::vector<std::string> volumeMaterials = splitString(line.at(6), ";");

    // The volume bounds shape
    std::string boundsShape = volumeBounds.at(0);

    // Split the bounds on the deliminater
    ACTS_DEBUG("Creating (" << volumeType << ") Blueprint node for volume "
                            << volumeName << " (id: " << volumeId
                            << ", shape: " << boundsShape << ")");

    // Create a table entry per defined volume
    TableEntry entry{volumeId,       volumeType,      volumeName,
                     volumeBounds,   volumeInternals, volumeBinnings,
                     volumeMaterials};

    // This will guarantee to have access to the building
    blueprintTableMap[volumeName] = entry;
  }

  // Now we can build the tree
  auto topEntry = blueprintTableMap.find(options.topEntry);
  if (topEntry == blueprintTableMap.end()) {
    throw std::invalid_argument("GeoModelBlueprintCreater: Top node '" +
                                options.topEntry +
                                "' not found in blueprint table");
  }

  return blueprint;
}

std::unique_ptr<Acts::Experimental::Blueprint::Node>
Acts::GeoModelBlueprintCreater::createNode(
    const TableEntry& entry,
    const std::vector<ActsScalar>& motherBounds) const {
  // Parse the bounds
  auto [boundsType, boundsValues, translation] =
      parseBounds(entry.bounds, motherBounds);

  Transform3 transform = Acts::Transform3::Identity();
  transform.translation() = translation;

  // Block for branch or container nodes that have children
  if (entry.type == "branch" || entry.type == "container" ||
      entry.type == "root") {
    std::vector<std::unique_ptr<Experimental::Blueprint::Node>> children;
    // Create the sub nodes
    for (const auto& childEntryStr : entry.internals) {
      auto childEntry = blueprintTableMap.find(childEntryStr);
      if (childEntry == blueprintTableMap.end()) {
        throw std::invalid_argument("GeoModelBlueprintCreater: Child node '" +
                                    childEntryStr + "' of '" + entry.name +
                                    "' found in blueprint table");
      }
      children.push_back(createNode(childEntry->second, boundsValues));
    }

    // Create the binnings
    std::vector<Acts::BinningValue> binnings;
    std::for_each(entry.binnings.begin(), entry.binnings.end(),
                  [&binnings](const std::string& b) {
                    binnings.push_back(toBinningValue(b));
                  });

    // Create the branch node
    return std::make_unique<Experimental::Blueprint::Node>(
        entry.name, transform, boundsType, boundsValues, binnings, children);

  } else if (entry.type == "leaf") {
    return std::make_unique<Experimental::Blueprint::Node>(
        entry.name, transform, boundsType, boundsValues, nullptr);
  } else {
    throw std::invalid_argument(
        "GeoModelBlueprintCreater: Unknown node type '" + entry.type + "'");
  }

  return nullptr;
}

std::tuple<Acts::VolumeBounds::BoundsType, std::vector<Acts::ActsScalar>,
           Acts::Vector3>
Acts::GeoModelBlueprintCreater::parseBounds(
    const std::vector<std::string>& boundsEntry,
    const std::vector<ActsScalar>& motherBounds) const {
  // Create the return values
  Acts::VolumeBounds::BoundsType boundsType;
  Acts::Vector3 translation{0., 0., 0.};
  std::vector<ActsScalar> boundsValues = {};

  // Switch on the bounds type
  if (boundsEntry[0u] == "cyl") {
    // Set the bounds type
    boundsType = Acts::VolumeBounds::BoundsType::eCylinder;
    // Capture the values
    std::vector<std::string> valuesEntry = splitString(boundsEntry[1u], ",");
    if (valuesEntry.size() < 4u) {
      throw std::invalid_argument(
          "GeoModelBlueprintCreater: Cylinder bounds entry has to have at "
          "least 4 entries (rmin, rmax, "
          "zmin, zmax)");
    }
    for (auto [iv, value] : Acts::enumerate(valuesEntry)) {
      ActsScalar val = 0.;
      if (value == "i") {
        if (iv >= motherBounds.size()) {
          throw std::invalid_argument(
              "GeoModelBlueprintCreater: Mother bounds do not have enough "
              "entries for value inheritance.");
        }
        val = motherBounds.at(iv);
      } else {
        val = std::stod(value);
      }
      boundsValues.push_back(val);
    }
    // Create the translation
    ActsScalar z = 0.5 * (boundsValues[2] + boundsValues[3]);
    translation = Acts::Vector3(0., 0., z);
  } else {
    throw std::invalid_argument(
        "GeoModelBlueprintCreater: Unknown bounds type, only 'cyl' is "
        "supported for the moment.");
  }

  return std::make_tuple(boundsType, boundsValues, translation);
}

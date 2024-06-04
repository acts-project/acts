// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelBlueprintCreater.hpp"

#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/detail/BlueprintHelper.hpp"
#include "Acts/Detector/interface/IGeometryIdGenerator.hpp"
#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include "Acts/Plugins/GeoModel/detail/GeoModelDbHelper.hpp"
#include "Acts/Plugins/GeoModel/detail/GeoModelExtentHelper.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/RangeXD.hpp"

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

}  // namespace

using namespace Acts::detail;

Acts::GeoModelBlueprintCreater::GeoModelBlueprintCreater(
    const Config& cfg, std::unique_ptr<const Logger> mlogger)
    : m_cfg(cfg), m_logger(std::move(mlogger)) {}

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

  // Prepare the map
  std::map<std::string, TableEntry> blueprintTableMap;

  Cache cache;
  // Prepare the KdtSurfaces if configured to do so
  //
  if (!m_cfg.detectorSurfaces.empty()) {
    std::array<BinningValue, 3u> kdtBinning = {binX, binY, binZ};
    if (m_cfg.kdtBinning.empty()) {
      throw std::invalid_argument(
          "GeoModelBlueprintCreater: At least one binning value for KDTree "
          "structure has to be given");
    } else if (m_cfg.kdtBinning.size() == 1u) {
      kdtBinning = {m_cfg.kdtBinning[0u], m_cfg.kdtBinning[0u],
                    m_cfg.kdtBinning[0u]};
    } else if (m_cfg.kdtBinning.size() == 2u) {
      kdtBinning = {m_cfg.kdtBinning[0u], m_cfg.kdtBinning[1u],
                    m_cfg.kdtBinning[1u]};
    } else if (m_cfg.kdtBinning.size() == 3u) {
      kdtBinning = {m_cfg.kdtBinning[0u], m_cfg.kdtBinning[1u],
                    m_cfg.kdtBinning[2u]};
    } else {
      throw std::invalid_argument(
          "GeoModelBlueprintCreater: Too many binning values for KDTree "
          "structure");
    }

    // Create the KdtSurfaces
    cache.kdtSurfaces = std::make_shared<Experimental::KdtSurfaces<3u>>(
        gctx, m_cfg.detectorSurfaces, kdtBinning);
  }

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
    // volume bounds of top volume might be overruled for top node
    std::vector<std::string> volumeBounds =
        (!options.topBoundsOverride.empty() && volumeName == options.topEntry)
            ? GeoModelDbHelper::splitString(options.topBoundsOverride, "|")
            : GeoModelDbHelper::splitString(line.at(3), "|");
    std::vector<std::string> volumeInternals =
        GeoModelDbHelper::splitString(line.at(4), ":");
    std::vector<std::string> volumeBinnings =
        GeoModelDbHelper::splitString(line.at(5), ";");
    std::vector<std::string> volumeMaterials =
        GeoModelDbHelper::splitString(line.at(6), ";");

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

  // Recursively create the nodes
  blueprint.name = topEntry->second.name;
  blueprint.topNode =
      createNode(cache, gctx, topEntry->second, blueprintTableMap, Extent());

  // Return the ready-to-use blueprint
  return blueprint;
}

std::unique_ptr<Acts::Experimental::Blueprint::Node>
Acts::GeoModelBlueprintCreater::createNode(
    Cache& cache, const GeometryContext& gctx, const TableEntry& entry,
    const std::map<std::string, TableEntry>& tableEntryMap,
    const Extent& motherExtent) const {
  ACTS_DEBUG("Build Blueprint node for '" << entry.name << "'.");

  // Peak into the volume entry to understand which one should be constraint
  // by the internals building
  auto internalContraints =
      detail::GeoModelExentHelper::readConstaints(entry.bounds, "i");

  // Create and return the container node with internal constrtins
  auto [internalsBuilder, internalExtent] = createInternalStructureBuilder(
      cache, gctx, entry, motherExtent, internalContraints);

  if (internalsBuilder != nullptr) {
    ACTS_VERBOSE("Internal building yielded extent "
                 << internalExtent.toString());
  }

  // Parse the bounds
  auto [boundsType, extent, boundValues, translation] =
      parseBounds(entry.bounds, motherExtent, internalExtent);

  ACTS_VERBOSE("Creating with extent " << extent.toString());

  Transform3 transform = Acts::Transform3::Identity();
  transform.translation() = translation;

  std::vector<std::string> entryTypeSplit =
      detail::GeoModelDbHelper::splitString(entry.type, "|");
  std::string entryType = entryTypeSplit[0u];

  // Block for branch or container nodes that have children
  if (entryType == "branch" || entryType == "container" ||
      entryType == "root") {
    std::vector<std::unique_ptr<Experimental::Blueprint::Node>> children;
    // Check for gap filling
    bool gapFilling = false;
    // Check if the entry has children
    if (entry.internals.size() < 2u) {
      throw std::invalid_argument(
          "GeoModelBlueprintCreater: Branch node '" + entry.name +
          "' has no children defined in blueprint table");
    }
    std::vector<std::string> childrenNames =
        GeoModelDbHelper::splitString(entry.internals[1u], ",");
    // Create the sub nodes and keep track of the raw values
    for (const auto& childName : childrenNames) {
      std::string fChildName = entry.name + std::string("/") + childName;
      if (childName == "*") {
        // Gap volume detected
        gapFilling = true;
        ACTS_VERBOSE("Gap volume detected, gap filling will be triggered.");
        continue;
      }
      // Check for child and build it
      auto childEntry = tableEntryMap.find(fChildName);
      if (childEntry == tableEntryMap.end()) {
        throw std::invalid_argument("GeoModelBlueprintCreater: Child node '" +
                                    childName + "' of '" + entry.name +
                                    "' NOT found in blueprint table");
      }
      children.push_back(
          createNode(cache, gctx, childEntry->second, tableEntryMap, extent));
    }

    // Create the binnings
    std::vector<Acts::BinningValue> binnings;
    std::for_each(entry.binnings.begin(), entry.binnings.end(),
                  [&binnings](const std::string& b) {
                    binnings.push_back(toBinningValue(b));
                  });

    // Complete the children
    auto node = std::make_unique<Experimental::Blueprint::Node>(
        entry.name, transform, boundsType, boundValues, binnings,
        std::move(children), extent);
    if (gapFilling) {
      // Find the first child that is not a gap
      Experimental::detail::BlueprintHelper::fillGaps(*node, true);
    }

    // Attach a geoID generator if configured
    if (entryTypeSplit.size() > 1u) {
      // Get the geometry ID from the string
      int geoID = std::stoi(entryTypeSplit[1u]);
      Experimental::GeometryIdGenerator::Config geoIDCfg;
      geoIDCfg.containerMode = true;
      geoIDCfg.containerId = geoID;
      geoIDCfg.resetSubCounters = true;
      // Make the container geoID generator
      node->geoIdGenerator = std::make_shared<Experimental::GeometryIdGenerator>(
          geoIDCfg, m_logger->clone(entry.name + "_GeometryIdGenerator"));
    }

    // Create the branch node
    return node;

  } else if (entryType == "leaf") {
    return std::make_unique<Experimental::Blueprint::Node>(
        entry.name, transform, boundsType, boundValues, internalsBuilder,
        extent);
  } else {
    throw std::invalid_argument(
        "GeoModelBlueprintCreater: Unknown node type '" + entry.type + "'");
  }

  return nullptr;
}

std::tuple<std::shared_ptr<const Acts::Experimental::IInternalStructureBuilder>,
           Acts::Extent>
Acts::GeoModelBlueprintCreater::createInternalStructureBuilder(
    Cache& cache, const GeometryContext& gctx, const TableEntry& entry,
    const Extent& externalExtent,
    const std::vector<BinningValue>& internalConstraints) const {
  // Check if the internals entry is empty
  if (entry.internals.empty()) {
    return std::make_tuple(nullptr, Extent());
  }

  // Build a layer structure
  if (entry.internals[0u] == "layer") {
    // Check if the internals entry is interpretable
    if (entry.internals.size() < 2u) {
      throw std::invalid_argument(
          "GeoModelBlueprintCreater: Internals entry not complete.");
    }

    // Internal split of the internals
    std::vector<std::string> internalsSplit =
        detail::GeoModelDbHelper::splitString(entry.internals[1u], "|");

    // Prepare an internal extent
    Extent internalExtent;
    if (internalsSplit[0u] == "kdt" && cache.kdtSurfaces != nullptr) {
      std::vector<std::string> internalsData = {internalsSplit[1u],
                                                internalsSplit[2u]};
      auto [boundsType, rangeExtent] =
          detail::GeoModelExentHelper::extentFromTable(internalsData,
                                                       externalExtent);

      ACTS_VERBOSE("Requested range: " << rangeExtent.toString());
      std::array<ActsScalar, 3u> mins = {};
      std::array<ActsScalar, 3u> maxs = {};

      // Fill what we have - follow the convention to fill up with the last
      for (std::size_t ibv = 0; ibv < 3u; ++ibv) {
        if (ibv < m_cfg.kdtBinning.size()) {
          BinningValue v = m_cfg.kdtBinning[ibv];
          mins[ibv] = rangeExtent.min(v);
          maxs[ibv] = rangeExtent.max(v);
          continue;
        }
        mins[ibv] = rangeExtent.min(m_cfg.kdtBinning.back());
        maxs[ibv] = rangeExtent.max(m_cfg.kdtBinning.back());
      }
      // Create the search range
      RangeXD<3u, ActsScalar> searchRange{mins, maxs};
      auto surfaces = cache.kdtSurfaces->surfaces(searchRange);
      // Loop over surfaces and create an internal extent
      for (auto& sf : surfaces) {
        auto sfExtent =
            sf->polyhedronRepresentation(gctx, m_cfg.nSegments).extent();
        internalExtent.extend(sfExtent, internalConstraints);
      }
      ACTS_VERBOSE("Found " << surfaces.size() << " surfaces in range "
                            << searchRange.toString());
      ACTS_VERBOSE("Internal extent: " << internalExtent.toString());

      // Create the layer structure builder
      Experimental::LayerStructureBuilder::Config lsbCfg;
      lsbCfg.surfacesProvider =
          std::make_shared<Experimental::LayerStructureBuilder::SurfacesHolder>(
              surfaces);
      lsbCfg.nMinimalSurfaces = surfaces.size() + 1u;

      return std::make_tuple(
          std::make_shared<Experimental::LayerStructureBuilder>(
              lsbCfg, m_logger->clone(entry.name + "_LayerStructureBuilder")),
          internalExtent);

    } else {
      throw std::invalid_argument(
          "GeoModelBlueprintCreater: Unknown layer internals type '" +
          entry.internals[1u] + "' / or now kdt surfaces provided.");
    }
  }
  return std::make_tuple(nullptr, Extent());
}

std::tuple<Acts::VolumeBounds::BoundsType, Acts::Extent,
           std::vector<Acts::ActsScalar>, Acts::Vector3>
Acts::GeoModelBlueprintCreater::parseBounds(
    const std::vector<std::string>& boundsEntry, const Extent& externalExtent,
    const Extent& internalExtent) const {
  // Create the return values
  Vector3 translation{0., 0., 0.};
  std::vector<ActsScalar> boundValues = {};
  auto [boundsType, extent] = detail::GeoModelExentHelper::extentFromTable(
      boundsEntry, externalExtent, internalExtent);

  // Switch on the bounds type
  if (boundsEntry[0u] == "cyl") {
    // Create the translation & bound values
    translation = Acts::Vector3(0., 0., extent.medium(binZ));
    boundValues = {extent.min(binR), extent.max(binR),
                   0.5 * extent.interval(binZ)};
  } else {
    throw std::invalid_argument(
        "GeoModelBlueprintCreater: Unknown bounds type, only 'cyl' is "
        "supported for the moment.");
  }

  return std::make_tuple(boundsType, extent, boundValues, translation);
}

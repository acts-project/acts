// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/DD4hepBlueprint.hpp"

#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/IndexedRootVolumeFinderBuilder.hpp"
#include "Acts/Plugins/DD4hep/DD4hepBinningHelpers.hpp"
#include "Acts/Plugins/DD4hep/DD4hepConversionHelpers.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <sstream>

Acts::Experimental::DD4hepBlueprint::DD4hepBlueprint(
    const Config& cfg, std::unique_ptr<const Logger> mlogger)
    : m_cfg(cfg), m_logger(std::move(mlogger)) {
  ACTS_DEBUG("UnitLength conversion factor (DD4hep -> Acts): " << unitLength);
}

std::unique_ptr<Acts::Experimental::Blueprint::Node>
Acts::Experimental::DD4hepBlueprint::create(
    Cache& cache, const dd4hep::DetElement& dd4hepElement) const {
  ACTS_DEBUG("Drawing a blueprint from the DD4hep element '"
             << dd4hepElement.name() << "'.");

  // Create the root node
  std::vector<ActsScalar> bValues = {0., 150., 1000.};
  std::vector<BinningValue> binning = {Acts::binR};
  auto root = std::make_unique<Acts::Experimental::Blueprint::Node>(
      dd4hepElement.name(), Acts::Transform3::Identity(),
      Acts::VolumeBounds::eCylinder, bValues, binning);

  // Recursively parse the tree
  recursiveParse(cache, *root, dd4hepElement);
  // Return the top node
  return root;
}

void Acts::Experimental::DD4hepBlueprint::recursiveParse(
    Cache& cache, Blueprint::Node& mother,
    const dd4hep::DetElement& dd4hepElement, unsigned int hiearchyLevel) const {
  // This will allow to skip empty hierarchy levels
  Blueprint::Node* current = &mother;
  unsigned int hierarchyAddOn = 0;

  std::string ofs(hiearchyLevel * 2u, ' ');

  // Node types
  std::vector<std::string> nodeTypes = {"acts_world", "acts_container",
                                        "acts_volume"};
  for (const auto& nType : nodeTypes) {
    // Check if it complies with the given definition
    bool ntt = getParamOr<bool>(nType, dd4hepElement, false);
    if (ntt) {
      ACTS_VERBOSE(ofs << "ACTS node '" << nType
                       << "' attached to dd4hep element '"
                       << dd4hepElement.name() << "',");
      // Extract the bounds type, values and binning
      auto [transform, bValueType, bValues, binning, auxExt] =
          extractExternals(dd4hepElement, nType);
      // Extract potential internal builders and tools
      auto [internalsBuilder, rootsFinderBuilder, geoIdGenerator, auxInt] =
          extractInternals(cache.dd4hepStore, dd4hepElement, nType);
      // Screen output of position and shape
      ACTS_VERBOSE(ofs << " - translation  : "
                       << toString(transform.translation()));
      ACTS_VERBOSE(ofs << " - bounds type  : " << bValueType);
      ACTS_VERBOSE(ofs << " - bound values : " << toString(bValues));
      // If it is not the world node, create a new one
      if (nType == "acts_world") {
        mother.transform = transform;
        mother.boundsType = bValueType;
        mother.boundaryValues = bValues;
        mother.binning = binning;

      } else if (nType == "acts_container") {
        // Creating the branch node
        auto branch = std::make_unique<Acts::Experimental::Blueprint::Node>(
            dd4hepElement.name(), transform, bValueType, bValues, binning);
        current = branch.get();
        mother.add(std::move(branch));

      } else if (nType == "acts_volume") {
        // Crreating a leaf node
        auto leaf = std::make_unique<Acts::Experimental::Blueprint::Node>(
            dd4hepElement.name(), transform, bValueType, bValues);
        current = leaf.get();
        mother.add(std::move(leaf));
      }
      // Current is set now appropriately, adding auxiliary information
      if (!auxExt.empty()) {
        ACTS_VERBOSE(ofs << " - " << auxExt);
        current->auxiliary.push_back(auxExt);
      }
      // Adding the internals builder - if present
      if (internalsBuilder != nullptr) {
        ACTS_VERBOSE(ofs << " - " << auxInt[0u]);
        current->internalsBuilder = internalsBuilder;
      }
      // Adding root finder builder - if present
      if (rootsFinderBuilder != nullptr) {
        ACTS_VERBOSE(ofs << " - " << auxInt[1u]);
        current->rootVolumeFinderBuilder = rootsFinderBuilder;
      }

      // Adding geo Id generator - if present
      if (geoIdGenerator != nullptr) {
        ACTS_VERBOSE(ofs << " - " << auxInt[2u]);
        current->geoIdGenerator = geoIdGenerator;
      }
    }
  }

  // Step down to the children - not possible for leaf nodes
  const dd4hep::DetElement::Children& children = dd4hepElement.children();
  if (!children.empty()) {
    ACTS_VERBOSE(ofs << "dd4hep element '" << dd4hepElement.name() << "' has "
                     << children.size() << " children.");
    for (auto& child : children) {
      dd4hep::DetElement dd4hepChild = child.second;
      recursiveParse(cache, *current, dd4hepChild,
                     hiearchyLevel + hierarchyAddOn);
    }
  }
}

std::tuple<Acts::Transform3, Acts::VolumeBounds::BoundsType,
           std::vector<Acts::ActsScalar>, std::vector<Acts::BinningValue>,
           std::string>
Acts::Experimental::DD4hepBlueprint::extractExternals(
    const dd4hep::DetElement& dd4hepElement,
    const std::string& baseName) const {
  std::string aux = "";
  // Get the bounds type
  auto bValueInt = Acts::getParamOr<int>(
      baseName + "_type", dd4hepElement,
      static_cast<int>(Acts::VolumeBounds::BoundsType::eOther));
  auto bValueType = static_cast<Acts::VolumeBounds::BoundsType>(bValueInt);
  // Get the bounds values
  auto bValues = Acts::extractSeries<Acts::ActsScalar>(
      dd4hepElement, baseName + "_bvalues", unitLength);
  // Get the binning values
  std::vector<Acts::BinningValue> bBinning = {};
  auto binningString =
      Acts::getParamOr<std::string>(baseName + "_binning", dd4hepElement, "");
  if (!binningString.empty()) {
    aux += "vol. binning : " + binningString;
    std::string del = ",";
    int end = binningString.find(del);
    end = (end == -1) ? binningString.length() : end;
    // Split and convert
    while (end != -1) {
      std::string b = binningString.substr(0, end);
      bBinning.push_back(Acts::stringToBinningValue(b));
      end = binningString.find(del);
    }
  }
  /// Get the transform
  auto transform = Acts::extractTransform(dd4hepElement, baseName, unitLength);
  // Return the tuple
  return std::make_tuple(transform, bValueType, bValues, bBinning, aux);
}

std::tuple<std::shared_ptr<const Acts::Experimental::IInternalStructureBuilder>,
           std::shared_ptr<const Acts::Experimental::IRootVolumeFinderBuilder>,
           std::shared_ptr<const Acts::Experimental::IGeometryIdGenerator>,
           std::array<std::string, 3u>>
Acts::Experimental::DD4hepBlueprint::extractInternals(
    Acts::DD4hepDetectorElement::Store& dd4hepStore,
    const dd4hep::DetElement& dd4hepElement,
    const std::string& baseName) const {
  // Return objects
  std::shared_ptr<const Acts::Experimental::IInternalStructureBuilder>
      internalsBuilder = nullptr;
  std::shared_ptr<const Acts::Experimental::IRootVolumeFinderBuilder>
      rootsFinderBuilder = nullptr;
  std::shared_ptr<const Acts::Experimental::IGeometryIdGenerator>
      geoIdGenerator = nullptr;
  /// Auxiliary information
  std::array<std::string, 3u> aux = {"", "", ""};

  // Check for internal structure builder
  auto internals =
      Acts::getParamOr<bool>(baseName + "_internals", dd4hepElement, false);
  if (internals) {
    auto internalsType = Acts::getParamOr<std::string>(
        baseName + "_internals_type", dd4hepElement, "");
    if (internalsType == "layer") {
      aux[0u] = "int. struct : layer";
      // Create a new layer builder
      DD4hepLayerStructure::Options lOptions;
      lOptions.name = dd4hepElement.name();
      internalsBuilder =
          m_cfg.layerStructure->builder(dd4hepStore, dd4hepElement, lOptions);
    }
  }

  // Check for root volume finder
  auto rootFinder = Acts::getParamOr<std::string>(
      baseName + "_root_volume_finder", dd4hepElement, "");
  if (rootFinder == "indexed") {
    aux[1u] = "root finder : indexed";
    std::vector<BinningValue> binning = {binZ, binR};
    rootsFinderBuilder =
        std::make_shared<Acts::Experimental::IndexedRootVolumeFinderBuilder>(
            binning);
  }

  // Check for geo Id generator
  auto geoIdGen =
      Acts::getParamOr<std::string>(baseName + "_geo_id", dd4hepElement, "");
  if (geoIdGen == "incremental") {
    aux[2u] = "geo_id gen. : incremental";
    Acts::Experimental::GeometryIdGenerator::Config geoIdCfg;
    geoIdGenerator =
        std::make_shared<Acts::Experimental::GeometryIdGenerator>(geoIdCfg);
  }

  return std::make_tuple(internalsBuilder, rootsFinderBuilder, geoIdGenerator,
                         aux);
}

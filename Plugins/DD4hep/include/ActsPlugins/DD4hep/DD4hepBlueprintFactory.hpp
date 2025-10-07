// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Detector/Blueprint.hpp"
#include "Acts/Detector/interface/IGeometryIdGenerator.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"
#include "Acts/Detector/interface/IRootVolumeFinderBuilder.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include "ActsPlugins/DD4hep/DD4hepLayerStructure.hpp"

#include <memory>
#include <optional>
#include <tuple>
#include <vector>

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetElement.h>

namespace ActsPlugins {

class DD4hepBlueprintFactory {
 public:
  /// @brief Nested config object
  struct Config {
    std::shared_ptr<ActsPlugins::DD4hepLayerStructure> layerStructure = nullptr;
    /// Configuration for DD4hep layer structure extraction

    /// The maximum number of portals to be checked for protal material
    unsigned int maxPortals = 8u;
  };

  /// @brief Nested cache object for the detector store
  struct Cache {
    /// Storage for DD4hep detector elements
    DD4hepDetectorElement::Store dd4hepStore;
  };

  /// @brief Constructor with arguments
  ///
  /// @param cfg the config object
  /// @param mlogger the logging instance
  explicit DD4hepBlueprintFactory(
      const Config& cfg,
      std::unique_ptr<const Acts::Logger> mlogger = Acts::getDefaultLogger(
          "DD4hepBlueprintFactory", Acts::Logging::INFO));

  /// @brief Create a blueprint from a DD4hep detector element
  ///
  /// @param cache the cache object for the detector store
  /// @param gctx the geometry context
  /// @param dd4hepElement the dd4hep detector element tree
  ///
  /// @return a new blueprint top node
  std::unique_ptr<Acts::Experimental::Gen2Blueprint::Node> create(
      Cache& cache, const Acts::GeometryContext& gctx,
      const dd4hep::DetElement& dd4hepElement) const;

 private:
  /// @brief auto-calculate the unit length conversion
  static constexpr double unitLength =
      Acts::UnitConstants::mm / dd4hep::millimeter;

  /// Configuration struct
  Config m_cfg;

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Private access to the logger
  const Acts::Logger& logger() const { return *m_logger; }

  /// @brief Recursive parsing of the detector element tree
  ///
  /// @param cache the cache object
  /// @param mother the mother node
  /// @param gctx the geometry context
  /// @param dd4hepElement the detector element at current level
  /// @param hiearchyLevel the current hierarchy level
  void recursiveParse(Cache& cache,
                      Acts::Experimental::Gen2Blueprint::Node& mother,
                      const Acts::GeometryContext& gctx,
                      const dd4hep::DetElement& dd4hepElement,
                      unsigned int hiearchyLevel = 0) const;

  /// @brief Extract the bounds type, values and binning from a DD4hep element
  ///
  /// @param gctx the geometry context
  /// @param dd4hepElement the DD4hep element
  /// @param baseName the base name of the acts type, e.g. "acts_volume", "acts_container", ...
  /// @param extOpt the optional extent as output from the internal parsing
  ///
  /// @return a tuple of the bounds type, values and binning, auxiliary data
  std::tuple<Acts::Transform3, Acts::VolumeBounds::BoundsType,
             std::vector<double>, std::vector<Acts::AxisDirection>, std::string>
  extractExternals(
      [[maybe_unused]] const Acts::GeometryContext& gctx,
      const dd4hep::DetElement& dd4hepElement, const std::string& baseName,
      const std::optional<Acts::Extent>& extOpt = std::nullopt) const;

  /// @brief Extract the internals from a DD4hep element
  ///
  /// This method parses the dd4hep element to return the internal structure
  /// builder, the root volume finder, and the geometry id generator.
  ///
  /// If any of the elements is not found, a nullptr is returned.
  ///
  /// @param dd4hepStore the store for the created dd4hep detector elements
  /// @param gctx the geometry context
  /// @param dd4hepElement the DD4hep element
  /// @param baseName the base name of the acts type, e.g. "acts_volume", "acts_container", ...
  ///
  /// @note The auxiliary information is returned as well for each of them
  ///
  /// @return a tuple of tools and auxiliary information
  std::tuple<
      std::shared_ptr<const Acts::Experimental::IInternalStructureBuilder>,
      std::shared_ptr<const Acts::Experimental::IRootVolumeFinderBuilder>,
      std::shared_ptr<const Acts::Experimental::IGeometryIdGenerator>,
      std::array<std::string, 3u>, std::optional<Acts::Extent>>
  extractInternals(DD4hepDetectorElement::Store& dd4hepStore,
                   const Acts::GeometryContext& gctx,
                   const dd4hep::DetElement& dd4hepElement,
                   const std::string& baseName) const;
};

}  // namespace ActsPlugins

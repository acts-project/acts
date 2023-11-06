// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Detector/Blueprint.hpp"
#include "Acts/Detector/interface/IGeometryIdGenerator.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"
#include "Acts/Detector/interface/IRootVolumeFinderBuilder.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Plugins/DD4hep/DD4hepLayerStructure.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <tuple>
#include <vector>

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetElement.h>

namespace Acts {
namespace Experimental {

class DD4hepBlueprint {
 public:
  /// @brief Nested config object
  struct Config {
    std::shared_ptr<Experimental::DD4hepLayerStructure> layerStructure =
        nullptr;
  };

  /// @brief Nested cache object for the detector store
  struct Cache {
    DD4hepDetectorElement::Store dd4hepStore;
  };

  /// @brief Constructor with arguments
  ///
  /// @param cfg the config object
  /// @param mlogger the logging instance
  DD4hepBlueprint(const Config& cfg,
                  std::unique_ptr<const Logger> mlogger =
                      getDefaultLogger("DD4hepBlueprint", Acts::Logging::INFO));

  /// @brief Create a blueprint from a DD4hep detector element
  ///
  /// @param cache the cache object for the detector store
  /// @param dd4hepElement the dd4hep detector element tree
  ///
  /// @return a new blueprint top node
  std::unique_ptr<Blueprint::Node> create(
      Cache& cache, const dd4hep::DetElement& dd4hepElement) const;

 private:
  /// @brief auto-calculate the unit length conversion
  static constexpr ActsScalar unitLength =
      Acts::UnitConstants::mm / dd4hep::millimeter;

  /// Configuration struct
  Config m_cfg;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to the logger
  const Logger& logger() const { return *m_logger; }

  /// @brief Recursive parsing of the detector element tree
  ///
  /// @param cache the cache object
  /// @param mother the mother node
  /// @param dd4hepElement the detector element at current level
  /// @param hiearchyLevel the current hierarchy level
  void recursiveParse(Cache& cache, Blueprint::Node& mother,
                      const dd4hep::DetElement& dd4hepElement,
                      unsigned int hiearchyLevel = 0) const;

  /// @brief Extract the bounds type, values and binning from a DD4hep element
  ///
  /// @param dd4hepElement the DD4hep element
  /// @param baseName the base name of the acts type, e.g. "acts_volume", "acts_container", ...
  ///
  /// @return a tuple of the bounds type, values and binning, auxiliary data
  std::tuple<Transform3, VolumeBounds::BoundsType, std::vector<ActsScalar>,
             std::vector<BinningValue>, std::string>
  extractExternals(const dd4hep::DetElement& dd4hepElement,
                   const std::string& baseName) const;

  /// @brief Extract the internals from a DD4hep element
  ///
  /// This method parses the dd4hep element to return the internal structure
  /// builder, the root volume finder, and the geometry id generator.
  ///
  /// If any of the elements is not found, a nullptr is returned.
  ///
  /// @param dd4hepStore the store for the created dd4hep detector elements
  /// @param dd4hepElement the DD4hep element
  /// @param baseName the base name of the acts type, e.g. "acts_volume", "acts_container", ...
  ///
  /// @note The auxiliary information is returned as well for each of them
  ///
  /// @return a tuple of tools and auxiliary information
  std::tuple<std::shared_ptr<const IInternalStructureBuilder>,
             std::shared_ptr<const IRootVolumeFinderBuilder>,
             std::shared_ptr<const IGeometryIdGenerator>,
             std::array<std::string, 3u>>
  extractInternals(DD4hepDetectorElement::Store& dd4hepStore,
                   const dd4hep::DetElement& dd4hepElement,
                   const std::string& baseName) const;
};

}  // namespace Experimental
}  // namespace Acts

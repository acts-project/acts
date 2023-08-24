// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/ProtoBinning.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <tuple>
#include <vector>

#include <DD4hep/DD4hepUnits.h>

class TGeoMatrix;

namespace dd4hep {
class DetElement;
}

namespace Acts {

class DD4hepDetectorElement;

/// A factory to convert DD4hep DetectorElements into sensitive
/// of passive surfaces
///
class DD4hepDetectorSurfaceFactory {
 public:
  /// Collect the sensitive surface & detector element
  using DD4hepSensitiveSurface =
      std::tuple<std::shared_ptr<DD4hepDetectorElement>,
                 std::shared_ptr<Surface>>;

  /// Collect the passive surfaces
  using DD4hepPassiveSurface = std::shared_ptr<Surface>;

  /// Collect passive surface proxies
  struct DD4hepPassiveSurfaceProxy {
    std::string placement = "";
  };

  /// Nested cache that records the conversion status
  struct Cache {
    /// The created detector elements - for the detector store
    std::vector<DD4hepSensitiveSurface> sensitiveSurfaces;
    /// The created non-const surfaces - for further processing,
    std::vector<DD4hepPassiveSurface> passiveSurfaces;
    /// The created Proxies for the support surfaces
    std::vector<DD4hepPassiveSurfaceProxy> passiveSurfaceProxies;
    /// matching and conversion statistics: surfaces
    std::size_t convertedSurfaces = 0;
    /// matching and conversion statistics: materials
    std::size_t convertedMaterials = 0;
    /// The collected binnings
    std::vector<Experimental::ProtoBinning> binnings = {};
    /// The collected supports
    std::vector<Experimental::LayerStructureBuilder::Support> supports = {};
  };

  /// The DD4hep detector element factory
  ///
  /// @param cfg the configuration struct
  /// @param mlogger a screen output logger
  DD4hepDetectorSurfaceFactory(
      std::unique_ptr<const Logger> mlogger = getDefaultLogger(
          "DD4hepDetectorSurfaceFactory", Acts::Logging::INFO));

  /// Construction method of the detector elements
  ///
  /// @param cache [in,out] into which the Elements are filled
  /// @param dd4hepElement the detector element representing the super structure
  ///
  /// @note this method will call the recursive construction
  void construct(Cache& cache, const dd4hep::DetElement& dd4hepElement);

 private:
  /// @brief  auto-calculate the unit length conversion
  static constexpr ActsScalar unitLength =
      Acts::UnitConstants::mm / dd4hep::millimeter;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  // The allows passive proxy positions
  std::vector<std::string> allowedPassiveProxies = {
      "inner", "outer", "negative", "positive", "representing"};

  /// Private access to the logger
  const Logger& logger() const { return *m_logger; }

  /// Construction method of the detector elements - recursive walk down
  ///
  /// @param cache [in,out] into which the Elements are filled
  /// @param dd4hepElement the detector element representing the super structure
  /// @param level the current level of the tree, used for log message output
  ///
  /// @note this method is called recursively
  void recursiveConstruct(Cache& cache, const dd4hep::DetElement& dd4hepElement,
                          int level);

  /// Method to convert a single sensitive detector element
  ///
  /// @param dd4hepElement the detector element
  ///
  /// @return a created detector element and surface
  DD4hepSensitiveSurface constructSensitiveElement(
      const dd4hep::DetElement& dd4hepElement) const;

  /// Method to convert a single sensitive detector element
  ///
  /// @param options the factory creation option
  /// @param dd4hepElement the detector element
  ///
  /// @return a created surface
  DD4hepPassiveSurface constructPassiveElement(
      const dd4hep::DetElement& dd4hepElement) const;
};

}  // namespace Acts

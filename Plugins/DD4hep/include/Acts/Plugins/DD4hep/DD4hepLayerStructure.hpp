// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/ProtoBinning.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorSurfaceFactory.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>

namespace dd4hep {
class DetElement;
}

namespace Acts {

namespace Experimental {

/// @brief This class allows to generate layer structure builders for dd4hep sub detectors
/// It performs an intermediate step by taking dd4hep::DetElemnent objects that
/// describe a detector structure and prepares the translation of the detector
/// element and eventual passive surfaces.
///
/// It would also build passive support structures if configured to do so.
///
class DD4hepLayerStructure {
 public:
  /// Constructor with DD4hepDetectorSurfaceFactory
  ///
  /// @param surfaceFactory the surfac factory which converts dd4hep::DetElement objects
  /// into their Acts coutnerparts
  /// @param logger is the screen output logger
  ///
  /// @note this needs to be provided
  DD4hepLayerStructure(
      std::shared_ptr<DD4hepDetectorSurfaceFactory> surfaceFactory,
      std::unique_ptr<const Logger> logger =
          getDefaultLogger("DD4hepLayerStructure", Acts::Logging::INFO));

  DD4hepLayerStructure() = delete;

  /// @brief Nested options struct
  ///
  /// If a binning description or a support cylinder
  /// description is chosen through options, it overwrites the corresponding
  /// description coming from DD4hep.
  struct Options {
    /// The name of the object
    std::string name = "";
    /// An out put log level
    Logging::Level logLevel = Logging::INFO;
    /// Approximation for the polyhedron binning nSegments
    unsigned int nSegments = 1u;
    /// Conversion options
    DD4hepDetectorSurfaceFactory::Options conversionOptions;
  };

  /// @brief This method generates a LayerStructure builder, which extends the
  /// IInternalStructureBuilder and can be used in the conjunction with a
  /// IExternalStructureBuilder to create `DetectorVolume` objects.
  ///
  /// It takes the detector element from DD4hep and some optional parameters
  ///
  /// @param dd4hepStore [in, out] the detector store for the built elements
  /// @param dd4hepElement the dd4hep detector element
  /// @param options containing the optional descriptions
  ///
  /// @return a LayerStructureBuilder
  std::shared_ptr<LayerStructureBuilder> builder(
      DD4hepDetectorElement::Store& dd4hepStore,
      const dd4hep::DetElement& dd4hepElement, const Options& options) const;

 private:
  /// The workorse: the surface factory
  std::shared_ptr<DD4hepDetectorSurfaceFactory> m_surfaceFactory = nullptr;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to the logger
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Experimental
}  // namespace Acts

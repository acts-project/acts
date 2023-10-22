// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Detector/VolumeStructureBuilder.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>

#include <DD4hep/DD4hepUnits.h>

namespace dd4hep {
class DetElement;
}

namespace Acts {

namespace Experimental {

/// @brief This class allows to generate volumes structure builders for dd4hep sub detectors,
/// together with an internal structure builder, this is sufficient to build the
/// new DetectorVolume objects.
///
class DD4hepVolumeStructure {
 public:
  /// Constructor with arguments
  /// @param mlogger is the screen output logger
  DD4hepVolumeStructure(std::unique_ptr<const Logger> mlogger =
                            getDefaultLogger("DD4hepVolumeStructure",
                                             Acts::Logging::INFO));

  DD4hepVolumeStructure() = delete;

  /// @brief nested options struct
  ///
  /// If a binning description or a support cylinder
  /// description is chosen through options, it overwrites the corresponding
  /// description coming from DD4hep.
  struct Options {
    /// The name of the object
    std::string name = "";
    /// An output log level for the builder
    Acts::Logging::Level logLevel = Acts::Logging::INFO;
  };

  /// @brief This method generates a VolumeStructureBuilder, which extends the
  /// IExternalStructureBuilder to create `DetectorVolume` objects.
  ///
  /// It takes the detector element from DD4hep and some optional parameters
  ///
  /// @param dd4hepElement the dd4hep detector element
  /// @param options containing the optional descriptions
  ///
  /// @return a VolumeStructureBuilder
  std::shared_ptr<VolumeStructureBuilder> builder(
      const dd4hep::DetElement& dd4hepElement, const Options& options) const;

 private:
  /// @brief  auto-calculate the unit length conversion
  static constexpr ActsScalar unitLength =
      Acts::UnitConstants::mm / dd4hep::millimeter;

  /// @brief  Method to recursively parse through the detector element
  /// and grep the volume structure information
  ///
  /// @param vsbConfig
  /// @param dd4hepElement
  ///
  /// @return a boolean break condition
  bool recursiveParse(VolumeStructureBuilder::Config& vsbConfig,
                      const dd4hep::DetElement& dd4hepElement) const;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to the logger
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Experimental
}  // namespace Acts

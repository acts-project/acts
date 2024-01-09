// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Digitization/DigitizationModule.hpp"

#include <memory>

#include "DD4hep/Segmentations.h"

namespace Acts {
class DigitizationModule;
}  // namespace Acts

/// In case several sensitive modules have the same segmentation it can and
/// should be shared between these modules to save memory and time.
/// In Acts the Acts::DigitizationModule is used to describe the geometric
/// digitization on a detector module. This Acts::DigitizationModule should be
/// shared amongst the modules with the same segmentation. In order to create it
/// there are currently two helper functions implemented
/// (ActsExamples::DD4hep::rectangleDigiModule(),ActsExamples::DD4hep::trapezoidalDigiModule)
/// which return the digitization module from DD4hep input. Afterwards an
/// ActsExtension from the same Acts::DigitizationModule can be created and
/// attached for all modules sharing the same segmentation.
///
/// Below you can find an example (in pseudo code) how to share the same
/// Acts::DigitizationModule
/// amongst modules (DetElements) which have the same segmentation in your
/// DD4hep detector constructor:
///
///  Create the Acts::DigitizationModule which should be shared amongst the
///  different modules using the global function with the dimensions of the
///  module and its DD4hep Segmentation. Where sensDet is the corresponding
///  DD4hep SensitiveDetector.
/// @code
/// auto digiModule = Acts::rectangularDigiModule(halflengthX,
///                                               halflnegthY,
///                                               thickness,
///                                               sensDet.readout().segmentation());
/// @endcode
/// Now loop over all modules which have the same segmentation,
/// create the Acts::ActsExtension from the digitization module
/// and attach the extension to the DD4hep::DetElement of the module (named
///   'moduleDetelement' here),
///
/// @code
/// for ('loop over modules') {
///   ...
///       Acts::ActsExtension* moduleExtension
///       = new Acts::ActsExtension(digiModule);
///   moduleDetElement.addExtension<Acts::ActsExtension>(moduleExtension);
/// }
/// @endcode
///
/// @param digiModule the Acts::DigitizationModule
/// @note in order to create the shared Acts::DigitizationModule from DD4hep
/// segmentation please use the global functions rectangleDigiModule() and
/// trapezoidalDigiModule().
///
/// If one wants to build the Acts Tracking Geometry with \a %DD4hep input
/// these
/// extension should be used during the construction of the \a %DD4hep
/// geometry
/// i.e. in the
/// \a %DD4hep detector constructors. First the ActsExtension configuration
/// object
/// should be created and then attached to the \a %DD4hep \a DetElement.
///
/// Example for a layer \a DetElement (\c layer_detElement) where also
/// parameters
/// for material mapping are handed over:
/// @code
///  Acts::ActsExtension::Config layConfig;
///  layConfig.isLayer               = true;
///  layConfig.axes                  = "XZy";
///  layConfig.materialBins1         = 50;
///  layConfig.materialBins2			= 100;
///  layConfig.layerMaterialPosition = Acts::LayerMaterialPos::inner
///  Acts::ActsExtension* layerExtension = new
///  Acts::ActsExtension(layConfig);
///  layer_detElement.addExtension<Acts::ActsExtension>(layerExtension);
///  @endcode
///
/// In case several sensitive detector modules have the same segmentation an
/// extension using the second constructor (with the segmentation as
/// parameter)
/// (or the function setSegmentation())
/// should be created once and then be attached to all the DetElements which
/// have that same segmentation. In this way only one
/// Acts::DigitizationModule
/// is
/// created and shared between all detector elements with the same
/// segmentation
/// which saves memory and time. If this extension is not set and the
/// DetElement
/// is sensitive and has a readout, a unique Acts::DigitizationModule will
/// be
/// created for this DetElement.
/// @endcode

namespace ActsExamples {
namespace DD4hep {

struct DD4hepDetectorHelper {
  /// Global method to build an Acts::DigitizationModule with rectangular
  /// segmentation.
  /// @note This function should be used in order to create the input
  /// needed for construction with
  /// Acts::ActsExtension(std::shared_ptr<const DigitizationModule>)
  /// @param halflengthX The half length in x of the detector module
  /// @param halflengthY The half length in y of the detector module
  /// @param thickness The thickness of the detector module
  /// @param segmentation the DD4hep segmentation
  static std::shared_ptr<const Acts::DigitizationModule> rectangleDigiModule(
      double halflengthX, double halflengthY, double thickness,
      const dd4hep::Segmentation& segmentation);

  /// Global method to build an Acts::DigitizationModule with trapezoidal
  /// segmentation.
  /// @note This function should be used in order to create the input
  /// needed for construction with
  /// Acts::ActsExtension(std::shared_ptr<const DigitizationModule>)
  /// @param minHalflengthX The half length in x of the detector module on the
  /// negative side of y
  /// @param maxHalflengthX The half length in x of the detector module on the
  /// positive side of y
  /// @param halflengthY The half length in y of the detector module
  /// @param thickness The thickness of the detector module
  /// @param segmentation the DD4hep segmentation
  static std::shared_ptr<const Acts::DigitizationModule> trapezoidalDigiModule(
      double minHalflengthX, double maxHalflengthX, double halflengthY,
      double thickness, const dd4hep::Segmentation& segmentation);
};

}  // end of namespace DD4hep
}  // end of namespace ActsExamples

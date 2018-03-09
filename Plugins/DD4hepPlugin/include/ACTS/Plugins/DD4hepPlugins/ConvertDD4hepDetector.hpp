// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_DD4HEPPLUGIN_CONVERTDD4HEPDETECTOR_H
#define ACTS_DD4HEPPLUGIN_CONVERTDD4HEPDETECTOR_H 1

#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Utilities/BinningType.hpp"
#include "ACTS/Utilities/Logger.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "DD4hep/DetElement.h"

namespace Acts {

/// @brief Global method which creates the TrackingGeometry from DD4hep input
///
/// This method returns a std::unique_ptr of the Acts::TrackingGeometry from the
/// World DD4hep \a DetElement.
/// @pre Before using this method make sure, that the preconditions described in
/// @ref DD4hepPlugins are met.
///
///
/// @param [in] worldDetElement the DD4hep DetElement of the world
/// @param [in] loggingLevel is the debug logging level of the conversion and
/// geometry building
/// @param [in] bTypePhi is how the sensitive surfaces (modules) should be
/// binned in
/// a layer in phi direction. Possible binningtypes:
/// 	- arbitrary   - of the sizes if the surfaces and the distance inbetween
/// 		vary. This mode finds out the bin boundaries by walking through the
/// 		surfaces.
/// 	- equidistant - if the sensitive surfaces are placed euqidistant
/// @note equidistant binningtype is recommended because it is faster not only
/// while building of the geometry but also for look up during the
/// extrapolation
/// @param [in] bTypeR is how the sensitive surfaces (modules) should be binned
/// in a
/// layer in r direction
/// 	- arbitrary   - if the sizes if the surfaces and the distance inbetween
/// 		them vary. This mode finds out the bin boundaries by walking through the
/// 		surfaces.
/// 	- equidistant - if the sensitive surfaces are placed euqidistant
/// @note equidistant binningtype is recommended because it is faster not only
/// while building of the geometry  but also for look up during the
/// extrapolation
/// @param [in] bTypeZ is how the sensitive surfaces (modules) should be binned
/// in a
/// layer in z direction
/// 	- arbitrary   - if the sizes if the surfaces and the distance inbetween
/// 		them vary. This mode finds out the bin boundaries by walking through the
/// 		surfaces.
/// 	- equidistant - if the sensitive surfaces are placed euqidistant
/// @note equidistant binningtype is recommended because it is faster not only
/// while building of the geometry  but also for look up during extrapolation
/// inner/outer
/// @param layerEnvelopeR the tolerance added to the geometrical extension in r
/// of the layers contained to build the volume envelope around
/// @param layerEnvelopeZ the tolerance added to the geometrical extension in z
/// of the layers contained to build the volume envelope around
/// @param buildDigitizationModules Flag indicating if the
/// Acts::DigitizationModule (needed for Acts geometric digitization) will be
/// build for every single sensitive DD4hep DetElement translating directly the
/// DD4hep Segmentation.
/// @attention Turning on this flag can be very time and memory consuming! If
/// different modules are sharing the same segmentation (which will be the case
/// most of the times) please use the
/// Acts::ActsExtension(std::shared_ptr<const DigitizationModule>) constructor.
/// More information on the usage can be found
/// in the description of the Acts::ActsExtension class.
/// @param defaultLayerThickness In case no surfaces (to be contained by the
/// layer) are handed over, the layer thickness will be set to this value
/// @note Layers containing surfaces per default are not allowed to be
///       attached to each other (navigation will fail at this point).
///       However, to allow material layers (not containing surfaces) to be
///       attached to each other, this default thickness is needed. In this
///       way, the layer will be thin (with space to the next layer), but
///       the material will have the'real' thickness.
/// @attention The default thickness should be set thin enough that no
///            touching or overlapping with the next layer can happen.
/// @exception std::logic_error if an error in the translation occurs
/// @return std::unique_ptr to the full Acts::TrackingGeometry
std::unique_ptr<const Acts::TrackingGeometry>
convertDD4hepDetector(dd4hep::DetElement worldDetElement,
                      Logging::Level     loggingLevel   = Logging::Level::INFO,
                      BinningType        bTypePhi       = equidistant,
                      BinningType        bTypeR         = equidistant,
                      BinningType        bTypeZ         = equidistant,
                      double             layerEnvelopeR = 1. * Acts::units::_mm,
                      double             layerEnvelopeZ = 1. * Acts::units::_mm,
                      bool               buildDigitizationModules = false,
                      double defaultLayerThickness = 10e-10 * Acts::units::_mm);

/// Method internally used by convertDD4hepDetector to collect all sub detectors
/// @param [in] detElement the dd4hep::DetElement of the volume of which the sub
/// detectors should be collected
/// @param [out] subdetectors the DD4hep::DetElements of the sub detectors
/// contained by detElement
void
collectSubDetectors(dd4hep::DetElement&              detElement,
                    std::vector<dd4hep::DetElement>& subdetectors);

/// Method internally used by convertDD4hepDetector to collect all volumes of a
/// compound detector
/// @param [in] detElement the dd4hep::DetElement of the volume of which the
/// compounds should be collected
/// @param [out] compounds the DD4hep::DetElements of the compounds contained by
/// detElement
void
collectCompounds(dd4hep::DetElement&              detElement,
                 std::vector<dd4hep::DetElement>& compounds);

/// Method internally used by convertDD4hepDetector
/// @param [in] detElement the dd4hep::DetElement of the volume of which the
/// layers should be collected
/// @param [out] layers the DD4hep::DetElements of the layers contained by
/// detElement
void
collectLayers(dd4hep::DetElement&              detElement,
              std::vector<dd4hep::DetElement>& layers);
}

#endif  // ACTS_DD4HEPPLUGIN_CONVERTDD4HEPDETECTOR_H

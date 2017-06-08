// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
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
#include "DD4hep/Detector.h"

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
/// @note this parameter only needs to be set if the volumes containing the
/// layers (e.g. barrel, endcap volumes) have no specific shape (assemblies)
/// @param layerEnvelopeZ the tolerance added to the geometrical extension in z
/// of the layers contained to build the volume envelope around
/// @note this parameter only needs to be set if the volumes containing the
/// layers (e.g. barrel, endcap volumes) have no specific shape (assemblies)
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
/// @exception std::logic_error if an error in the translation occurs
/// @return std::unique_ptr to the full Acts::TrackingGeometry
std::unique_ptr<const Acts::TrackingGeometry>
convertDD4hepDetector(DD4hep::Geometry::DetElement worldDetElement,
                      Logging::Level loggingLevel   = Logging::Level::INFO,
                      BinningType    bTypePhi       = equidistant,
                      BinningType    bTypeR         = equidistant,
                      BinningType    bTypeZ         = equidistant,
                      double         layerEnvelopeR = 0.,
                      double         layerEnvelopeZ = 0.,
                      bool           buildDigitizationModules = false);

/// Method internally used by convertDD4hepDetector
/// @param [in] detElement the DD4hep::DetElement of the volume of which the
/// layers should be collected
/// @param [out] layers the DD4hep::DetElements of the layers contained by
/// detElement
void
collectLayers(DD4hep::Geometry::DetElement&              detElement,
              std::vector<DD4hep::Geometry::DetElement>& layers);
}

#endif  // ACTS_DD4HEPPLUGIN_CONVERTDD4HEPDETECTOR_H

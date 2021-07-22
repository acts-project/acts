// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <DD4hep/DetElement.h>
#include <DD4hep/DetFactoryHelper.h>

namespace Acts {

class ProtoSurfaceMaterial;
class Layer;

/// Helper method to translate DD4hep material to Acts::ISurfaceMaterial
///
/// This is used to assign proto material to Cylinder Layers
///
/// @param detElement the DD4hep detector element (source)
/// @param cylinderLayer the Layer to be decorated (target)
/// @param loggingLevel is the output level for the conversion
///
/// @return a map of the identification string and a surface material
void addCylinderLayerProtoMaterial(
    dd4hep::DetElement detElement, Layer& cylinderLayer,
    Logging::Level loggingLevel = Logging::Level::INFO);

/// Helper method to translate DD4hep material to Acts::ISurfaceMaterial
///
/// Thisis used to assign proto material to Disc Layers
///
/// @param detElement the DD4hep detector element (source)
/// @param discLayer the Layer to be decorated (target)
/// @param loggingLevel is the output level for the conversion
///
/// @return a map of the identification string and a surface material
void addDiscLayerProtoMaterial(
    dd4hep::DetElement detElement, Layer& discLayer,
    Logging::Level loggingLevel = Logging::Level::INFO);

/// Helper method to be called for Cylinder and Disc Proto material
///
/// For both, cylinder and disc, the closed binning value is "binPhi"
///
/// @param aExtension the ActsExtension for the binning parameters
/// @param layer the Layer to assign the proto material (target)
/// @param binning the Binning prescription for the ActsExtension
void addLayerProtoMaterial(
    const ActsExtension& actsExtension, Layer& layer,
    const std::vector<std::pair<const std::string, Acts::BinningOption> >&
        binning);

/// Helper method to create proto material - to be called from the
/// addProto(...) methods
///
/// @param actsExtension the ActExtension to be checked
/// @param valueTag the xml tag for to ActsExtension to be parsed
/// @param firstBinning string lookup for first bin
/// @param binning the Binning prescription for the ActsExtension
std::shared_ptr<Acts::ProtoSurfaceMaterial> createProtoMaterial(
    const ActsExtension& actsExtension, const std::string& valueTag,
    const std::vector<std::pair<const std::string, Acts::BinningOption> >&
        binning);

/// Helper method that decorates an ActsExtension with proto material
/// description for boundaries
/// - it assigns bins for inner / representing / outer

/// Helper method that decorates an ActsExtension with proto material
/// description,
/// - it assigns bins for inner / representing / outer
///
/// @param x_material the material tag to be inspected
/// @param actsExtension the extension that is augmented
/// @param baseTag the xml tag to be checked
void xmlToProtoSurfaceMaterial(const xml_comp_t& x_material,
                               ActsExtension& actsExtension,
                               const std::string& baseTag);

}  // namespace Acts

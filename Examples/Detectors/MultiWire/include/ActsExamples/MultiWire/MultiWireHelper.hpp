// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/DetectorVolumeBuilder.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/VolumeStructureBuilder.hpp"
#include "Acts/Detector/interface/IExternalStructureBuilder.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <algorithm>
#include <iostream>
#include <vector>

namespace ActsExamples {

namespace MultiWireHelper {

/// Global method to get straw type surfaces from a gdml file by providing only
/// the names of sensitive and passive elements
/// @param sensitiveNames The name of sensitive elements for the Name Selector
/// @param passiveNames The name of passive elements for the Name Selector

std::vector<std::shared_ptr<Acts::Surface>> getStrawSurfaces(
    std::vector<std::string> sensitiveNames,
    std::vector<std::string> passiveNames);

/// Global method to create a 2-dimensional layer structure builder binning
/// providing the surfaces and the bounds for the grid assigned to the surfaces
/// @param surfaces The surfaces of the grid
/// @param multiWireBounds The bounds of the grid
std::vector<Acts::Experimental::LayerStructureBuilder::Binning> layerBinning(
    std::vector<std::shared_ptr<Acts::Surface>> surfaces,
    std::array<std::pair<float, float>, 3> multiWireBounds);

/// Global method to get the bounds of a multilayer with wires
/// @param surfaces The surfaces (wires) of the multilayer

std::array<std::pair<float, float>, 3> getMultiWireBounds(
    std::vector<std::shared_ptr<Acts::Surface>> surfaces);

/// Global method to create an internal strucutre builder for a multilayer using
/// the Layer Structure Builder extension
/// @param mlConfig The configure of the Layer Structure Builder
std::shared_ptr<Acts::Experimental::LayerStructureBuilder> internalLayerBuilder(
    Acts::Experimental::LayerStructureBuilder::Config mlConfig);

/// Global method to create an external structure builder for a multilayer using
/// the Volume Structure Builder extension
/// @param vsConfig The configure of the Volume Structure Builder

std::shared_ptr<Acts::Experimental::VolumeStructureBuilder>
externalVolumeBuilder(
    Acts::Experimental::VolumeStructureBuilder::Config vsConfig);

}  // namespace MultiWireHelper
}  // namespace ActsExamples

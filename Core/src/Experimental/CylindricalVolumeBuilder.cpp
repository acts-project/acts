// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Experimental/CylindricalVolumeBuilder.hpp"

std::vector<std::shared_ptr<Acts::DetectorVolume>>
Acts::CylindricalVolumeBuilder::volumesInZ(
    const std::vector<LayerBlueprint>& protoLayers,
    const std::string& name) {
  return {};
}

std::vector<std::shared_ptr<Acts::DetectorVolume>>
Acts::CylindricalVolumeBuilder::volumesInR(
    const std::vector<LayerBlueprint>& protoLayers,
    const std::string& name) {
  return {};
}
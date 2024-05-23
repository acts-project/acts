// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <GeoModelRead/ReadGeoModel.h>

class GeoVPhysVol;

namespace Acts {

struct GeoModelTree {
  std::shared_ptr<GeoModelIO::ReadGeoModel> geoReader = nullptr;
  GeoVPhysVol* worldVolume = nullptr;
  std::string worldVolumeName = "World";
};
}  // namespace Acts

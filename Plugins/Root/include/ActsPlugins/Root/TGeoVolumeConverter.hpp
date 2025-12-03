// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

class TGeoShape;
class TGeoMatrix;

namespace Acts {
class Volume;
}

namespace ActsPlugins {
struct TGeoVolumeConverter {
  TGeoVolumeConverter() = delete;

  static std::unique_ptr<Acts::Volume> cylinderVolume(
      const TGeoShape& tgShape, const TGeoMatrix& tgTransform,
      double lengthScale = 10.);
};

}  // namespace ActsPlugins

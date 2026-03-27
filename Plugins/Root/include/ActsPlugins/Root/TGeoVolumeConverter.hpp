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
/// Helper struct for converting supported ROOT `TGeoShape` instances into ACTS
/// volumes.
struct TGeoVolumeConverter {
  /// Utility type; not instantiable.
  TGeoVolumeConverter() = delete;

  /// Convert a cylindrical ROOT shape into an ACTS volume.
  /// @param tgShape ROOT shape to convert.
  /// @param tgTransform ROOT transform describing the placed volume.
  /// @param lengthScale Unit scale applied when converting ROOT lengths to
  ///        ACTS units.
  /// @return Newly created ACTS volume.
  static std::unique_ptr<Acts::Volume> cylinderVolume(
      const TGeoShape& tgShape, const TGeoMatrix& tgTransform,
      double lengthScale = 10.);
};

}  // namespace ActsPlugins

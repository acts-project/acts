// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/AbstractVolume.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Utilities/BoundingBox.hpp"

#include <array>
#include <vector>

namespace Acts {
namespace Test {

struct CubicBVHTrackingGeometry {
  using Box = Acts::Volume::BoundingBox;

  /// Default constructor for the Cubic Bounding Volume Hierarchy tracking
  /// geometry
  ///
  /// @param n number of boxes in each direction
  /// @param hl Range of the volume
  /// @param octd maximum depth
  CubicBVHTrackingGeometry(size_t n = 29, double hl = 1000, size_t octd = 5) {
    Box::Size size(Acts::Vector3(2, 2, 2));

    std::shared_ptr<CuboidVolumeBounds> vbds =
        std::make_shared<CuboidVolumeBounds>(10, 10, 10);

    double min = -hl;
    double max = hl;

    double step = (max - min) / double(n);
    std::vector<std::unique_ptr<const Volume>> boxVolumes;
    std::vector<std::unique_ptr<Box>> boxStore;
    boxStore.reserve((n + 1) * (n + 1) * (n + 1));

    std::vector<Box*> boxes;
    boxes.reserve(boxStore.size());

    for (size_t i = 0; i <= n; i++) {
      for (size_t j = 0; j <= n; j++) {
        for (size_t k = 0; k <= n; k++) {
          Vector3 pos(min + i * step, min + j * step, min + k * step);

          auto trf = Transform3(Translation3(pos));
          auto vol = std::make_unique<AbstractVolume>(trf, vbds);

          boxVolumes.push_back(std::move(vol));
          boxStore.push_back(
              std::make_unique<Box>(boxVolumes.back()->boundingBox()));
          boxes.push_back(boxStore.back().get());
        }
      }
    }

    Box* top = make_octree(boxStore, boxes, octd);

    // create trackingvolume
    // will own the volumes, so make non-owning copy first
    volumes.reserve(boxVolumes.size());
    for (auto& vol : boxVolumes) {
      volumes.push_back(vol.get());
    }

    // box like overall shape
    auto tvBounds =
        std::make_shared<CuboidVolumeBounds>(hl * 1.1, hl * 1.1, hl * 1.1);

    auto tv = TrackingVolume::create(Transform3::Identity(), tvBounds,
                                     std::move(boxStore), std::move(boxVolumes),
                                     top, nullptr, "TheVolume");

    trackingGeometry = std::make_shared<TrackingGeometry>(tv);
  }
  std::vector<const Volume*> volumes;
  std::shared_ptr<TrackingGeometry> trackingGeometry;
};

}  // namespace Test
}  // namespace Acts

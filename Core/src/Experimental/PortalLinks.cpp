// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Experimental/VolumeLink.hpp"

Acts::VolumeLink::VolumeLink(std::shared_ptr<DetectorVolume> volume,
                             std::shared_ptr<SurfaceArray> surfaceArray)
    : m_volumes({std::move(volume)}),
      m_surfaceArrays({std::move(surfaceArray)}),
      m_link(DefaultLink{}) {}

Acts::VolumeLink::VolumeLink(std::vector<DetectorVolumePtr> volumes,
                             std::vector<SurfaceArrayPtr> surfaceArrays,
                             Link&& link)
    : m_volumes(volumes),
      m_surfaceArrays(surfaceArrays),
      m_link(std::move(link)) {}

std::tuple<const Acts::DetectorVolume*, const Acts::SurfaceArray*>
Acts::VolumeLink::resolve(const GeometryContext& gctx, const Vector3& pos,
                          const Vector3& mom, NavigationDirection nDir) const {
  // Get the signed direction & get the corresponding bin
  const unsigned int bin = m_link(pos, dir);
  // Get the associated volumes
  const DetectorVolume* dVolume = nullptr;
  if (not m_volumes.empty()) {
    dVolume = m_volumes[bin].get();
  }
  // Get the associated surface array
  const SurfaceArray* sArray = nullptr;
  if (not m_surfaceArrays.empty()) {
    sArray = m_surfaceArrays[bin].get();
  }
  return {dVolume, sArray};
}

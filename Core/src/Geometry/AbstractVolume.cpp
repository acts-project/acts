// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/AbstractVolume.hpp"

#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <utility>

Acts::AbstractVolume::AbstractVolume(
    const Transform3& transform, std::shared_ptr<const VolumeBounds> volbounds)
    : Volume(transform, std::move(volbounds)) {
  createBoundarySurfaces();
}

const std::vector<Acts::BoundarySurfacePtr>&
Acts::AbstractVolume::boundarySurfaces() const {
  return m_boundarySurfaces;
}

void Acts::AbstractVolume::createBoundarySurfaces() {
  using Boundary = BoundarySurfaceT<AbstractVolume>;

  // Transform Surfaces To BoundarySurfaces
  auto orientedSurfaces = Volume::volumeBounds().orientedSurfaces(m_transform);

  m_boundarySurfaces.reserve(orientedSurfaces.size());
  for (auto& osf : orientedSurfaces) {
    AbstractVolume* opposite = nullptr;
    AbstractVolume* along = nullptr;
    if (osf.second == Direction::Negative) {
      opposite = this;
    } else {
      along = this;
    }
    m_boundarySurfaces.push_back(std::make_shared<const Boundary>(
        std::move(osf.first), opposite, along));
  }
}

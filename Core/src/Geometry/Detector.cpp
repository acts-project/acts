// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Detector.hpp"

#include "Acts/Geometry/NavigationState.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"

Acts::Experimental::Detector::Detector(
    const std::string& name,
    const std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    ManagedDetectorVolumeFinder&& volumeFinder)
    : m_name(name), m_volumeFinder(std::move(volumeFinder)) {
  if (volumes.empty()) {
    throw std::invalid_argument("Detector: no volumes were given.");
  }
  if (not volumeFinder.delegate.connected()) {
    throw std::invalid_argument(
        "Detector: volume finder delegate is not connected.");
  }
  // Fill and make unique
  std::vector<std::shared_ptr<DetectorVolume>> uniqueVolumes = volumes;
  for (auto v : volumes) {
    for (auto vv : v->volumePtrs()) {
      uniqueVolumes.push_back(vv);
    }
  }
  // Only keep the unique ones & fill the volume store
  auto last = std::unique(uniqueVolumes.begin(), uniqueVolumes.end());
  uniqueVolumes.erase(last, uniqueVolumes.end());
  m_volumes = DetectorVolume::ObjectStore<std::shared_ptr<DetectorVolume>>(
      uniqueVolumes);

  // Check for unique names and fill the volume name / index map
  for (auto [iv, v] : enumerate(m_volumes.external)) {
    const std::string vName = v->name();
    if (m_volumeNameIndex.find(vName) != m_volumeNameIndex.end()) {
      throw std::invalid_argument("Detector: duplicate volume name " + vName +
                                  " detected.");
    }
    m_volumeNameIndex[vName] = iv;
  }
}

std::shared_ptr<Acts::Experimental::Detector>
Acts::Experimental::Detector::getSharedPtr() {
  return shared_from_this();
}

std::shared_ptr<const Acts::Experimental::Detector>
Acts::Experimental::Detector::getSharedPtr() const {
  return shared_from_this();
}

const Acts::Experimental::DetectorVolume*
Acts::Experimental::Detector::findVolume(const GeometryContext& gctx,
                                         const Vector3& position) const {
  return m_volumeFinder.delegate(gctx, *this, position);
}

const Acts::Experimental::DetectorVolume*
Acts::Experimental::Detector::findVolume(const std::string& name) const {
  auto vCandidate = m_volumeNameIndex.find(name);
  if (vCandidate != m_volumeNameIndex.end()) {
    return m_volumes.external[vCandidate->second];
  }
  return nullptr;
}

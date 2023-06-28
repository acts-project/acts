// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Detector.hpp"

#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <iterator>
#include <stdexcept>
#include <utility>

Acts::Experimental::Detector::Detector(
    const std::string& name,
    std::vector<std::shared_ptr<DetectorVolume>> rootVolumes,
    DetectorVolumeUpdator&& detectorVolumeUpdator)
    : m_name(name),
      m_rootVolumes(std::move(rootVolumes)),
      m_detectorVolumeUpdator(std::move(detectorVolumeUpdator)) {
  if (m_rootVolumes.internal.empty()) {
    throw std::invalid_argument("Detector: no volume were given.");
  }
  if (not m_detectorVolumeUpdator.connected()) {
    throw std::invalid_argument(
        "Detector: volume finder delegate is not connected.");
  }

  // Fill volumes
  auto collectVolumes = [&]() {
    std::vector<std::shared_ptr<DetectorVolume>> volumes;
    auto recurse = [&volumes](const std::shared_ptr<DetectorVolume>& volume,
                              auto& callback) -> void {
      volumes.push_back(volume);
      for (const auto& v : volume->volumePtrs()) {
        callback(v, callback);
      }
    };
    for (const auto& root : m_rootVolumes.internal) {
      recurse(root, recurse);
    }
    return volumes;
  };
  m_volumes = DetectorVolume::ObjectStore<std::shared_ptr<DetectorVolume>>(
      collectVolumes());

  // Check for unique names and fill the volume name / index map
  for (auto [iv, v] : enumerate(m_volumes.internal)) {
    // Assign this detector
    v->assignDetector(*this);
    // Close the portals
    v->closePortals();
    // Store the name
    const std::string vName = v->name();
    if (m_volumeNameIndex.find(vName) != m_volumeNameIndex.end()) {
      throw std::invalid_argument("Detector: duplicate volume name " + vName +
                                  " detected.");
    }
    m_volumeNameIndex[vName] = iv;
  }
}

std::shared_ptr<Acts::Experimental::Detector>
Acts::Experimental::Detector::makeShared(
    const std::string& name,
    std::vector<std::shared_ptr<DetectorVolume>> rootVolumes,
    DetectorVolumeUpdator&& detectorVolumeUpdator) {
  return std::shared_ptr<Detector>(new Detector(
      name, std::move(rootVolumes), std::move(detectorVolumeUpdator)));
}

std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>&
Acts::Experimental::Detector::rootVolumePtrs() {
  return m_rootVolumes.internal;
}

const std::vector<const Acts::Experimental::DetectorVolume*>&
Acts::Experimental::Detector::rootVolumes() const {
  return m_rootVolumes.external;
}

std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>&
Acts::Experimental::Detector::volumePtrs() {
  return m_volumes.internal;
}

const std::vector<const Acts::Experimental::DetectorVolume*>&
Acts::Experimental::Detector::volumes() const {
  return m_volumes.external;
}

void Acts::Experimental::Detector::updateDetectorVolumeFinder(
    DetectorVolumeUpdator&& detectorVolumeUpdator) {
  m_detectorVolumeUpdator = std::move(detectorVolumeUpdator);
}

const Acts::Experimental::DetectorVolumeUpdator&
Acts::Experimental::Detector::detectorVolumeFinder() const {
  return m_detectorVolumeUpdator;
}

const std::string& Acts::Experimental::Detector::name() const {
  return m_name;
}

std::shared_ptr<Acts::Experimental::Detector>
Acts::Experimental::Detector::getSharedPtr() {
  return shared_from_this();
}

std::shared_ptr<const Acts::Experimental::Detector>
Acts::Experimental::Detector::getSharedPtr() const {
  return shared_from_this();
}

void Acts::Experimental::Detector::updateDetectorVolume(
    const GeometryContext& gctx, NavigationState& nState) const {
  m_detectorVolumeUpdator(gctx, nState);
}

const Acts::Experimental::DetectorVolume*
Acts::Experimental::Detector::findDetectorVolume(
    const GeometryContext& gctx, const Vector3& position) const {
  NavigationState nState;
  nState.currentDetector = this;
  nState.position = position;
  m_detectorVolumeUpdator(gctx, nState);
  return nState.currentVolume;
}

const Acts::Experimental::DetectorVolume*
Acts::Experimental::Detector::findDetectorVolume(
    const std::string& name) const {
  auto vCandidate = m_volumeNameIndex.find(name);
  if (vCandidate != m_volumeNameIndex.end()) {
    return m_volumes.external[vCandidate->second];
  }
  return nullptr;
}

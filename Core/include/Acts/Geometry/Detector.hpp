// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationDelegates.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace Acts {

namespace Experimental {

class Detector : public std::enable_shared_from_this<Detector> {
 protected:
  /// Create a detector from volumes
  ///
  /// @param name the detecor name
  /// @param volumes the objets contained by this detector
  /// @param volumeFinder is a Delegate to find the assocaited volume
  ///
  /// @note will throw an exception if volumes vector is empty
  /// @note will throw an exception if duplicate volume names exist
  /// @note will throw an exception if the delegate is not connected
  Detector(const std::string& name,
           const std::vector<std::shared_ptr<DetectorVolume>>& volumes,
           DetectorVolumeUpdator&& volumeFinder) noexcept(false);

 public:
  /// Factory for producing memory managed instances of Detector.
  /// Will forward all parameters and will attempt to find a suitable
  /// constructor.
  ///
  /// @tparam Args the arguments that will be forwarded
  template <typename... Args>
  static std::shared_ptr<Detector> makeShared(Args&&... args) {
    return std::shared_ptr<Detector>(new Detector(std::forward<Args>(args)...));
  }

  /// Retrieve a @c std::shared_ptr for this surface (non-const version)
  ///
  /// @note Will error if this was not created through the @c makeShared factory
  ///       since it needs access to the original reference. In C++14 this is
  ///       undefined behavior (but most likely implemented as a @c bad_weak_ptr
  ///       exception), in C++17 it is defined as that exception.
  /// @note Only call this if you need shared ownership of this object.
  ///
  /// @return The shared pointer
  std::shared_ptr<Detector> getSharedPtr();

  /// Retrieve a @c std::shared_ptr for this surface (const version)
  ///
  /// @note Will error if this was not created through the @c makeShared factory
  ///       since it needs access to the original reference. In C++14 this is
  ///       undefined behavior, but most likely implemented as a @c bad_weak_ptr
  ///       exception, in C++17 it is defined as that exception.
  /// @note Only call this if you need shared ownership of this object.
  ///
  /// @return The shared pointer
  std::shared_ptr<const Detector> getSharedPtr() const;

  /// Non-const access to the volumes
  ///
  /// @return the volumes shared pointer store
  std::vector<std::shared_ptr<DetectorVolume>>& volumePtrs();

  /// Const access to sub volumes
  ///
  /// @return a vector to const DetectorVolume raw pointers
  const std::vector<const DetectorVolume*>& volumes() const;

  /// Update the current volume of a given navigation state
  ///
  /// @param gctx is the Geometry context of the call
  /// @param nState [in, out] is the navigation state
  ///
  void updateDetectorVolume(const GeometryContext& gctx,
                            NavigationState& nState) const;

  /// Find a volume from a position
  ///
  /// @param gctx is the Geometry context of the call
  /// @param position is the position of the call
  ///
  /// @note this creates internally a NavigationState object
  ///
  /// @return the volume pointer or nullptr (if outside)
  const DetectorVolume* findDetectorVolume(const GeometryContext& gctx,
                                           const Vector3& position) const;

  /// Find a volume by name
  ///
  /// @param name with which the volume is searched for
  ///
  /// @return the volume pointer or nullptr (if not found)
  const DetectorVolume* findDetectorVolume(const std::string& name) const;

  /// Update the volume finder
  ///
  /// @param mVolumeFinder the new volume finder
  void updateDetectorVolumeFinder(DetectorVolumeUpdator&& mVolumeFinder);

  /// Const access to the volume finder
  const DetectorVolumeUpdator& detectorVolumeFinder() const;

  /// Return the name of the detector
  const std::string& name() const;

 private:
  /// Name of the detector
  std::string m_name = "Unnamed";

  /// Volume store (internal/external)
  DetectorVolume::ObjectStore<std::shared_ptr<DetectorVolume>> m_volumes;

  /// A volume finder delegate
  DetectorVolumeUpdator m_volumeFinder;

  /// Name/index map to find volumes by name and detect duplicates
  std::unordered_map<std::string, size_t> m_volumeNameIndex;
};

inline std::vector<std::shared_ptr<DetectorVolume>>& Detector::volumePtrs() {
  return m_volumes.internal;
}

inline const std::vector<const DetectorVolume*>& Detector::volumes() const {
  return m_volumes.external;
}

inline void Detector::updateDetectorVolumeFinder(
    DetectorVolumeUpdator&& mVolumeFinder) {
  m_volumeFinder = std::move(mVolumeFinder);
}

inline const DetectorVolumeUpdator& Detector::detectorVolumeFinder() const {
  return m_volumeFinder;
}

inline const std::string& Detector::name() const {
  return m_name;
}

}  // namespace Experimental

}  // namespace Acts

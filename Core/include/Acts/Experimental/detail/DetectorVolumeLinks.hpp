// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <exception>
#include <memory>
#include <tuple>
#include <vector>

namespace Acts {
namespace Experimental {

class DetectorVolume;

namespace detail {

/// @brief a null volume link - explicitely
///
/// @note the method parameters are ignored
///
/// @return the registred volume
inline static const DetectorVolume* nullVolumeLink(
    const GeometryContext& /*ignored*/, const Vector3& /*ignored*/,
    const Vector3& /*ignored*/) {
  return nullptr;
}

/// @brief  The implementation of a single link,
/// The given volume is returned and geometry context,
/// position and direction ignored.
class SingleLinkImpl : public INavigationDelegate {
 public:
  /// The single volume to point to
  const DetectorVolume* dVolume = nullptr;
  /// Convencience constructor
  ///
  /// @param detectorVolume is the volume you point to
  SingleLinkImpl(const DetectorVolume& detectorVolume)
      : dVolume(&detectorVolume) {}

  /// @brief  the function call to be connected with the Delegete
  ///
  /// @note the method parameters are ignored
  ///
  /// @return the registred volume
  const DetectorVolume* targetVolume(const GeometryContext& /*ignored*/,
                                     const Vector3& /*ignored*/,
                                     const Vector3& /*ignored*/) const {
    return dVolume;
  }
};

/// @brief  The implementation of a multi link relationship
/// based on 1D information
///
/// The given volume is returned and geometry context,
/// position and direction ignored.
///
/// @note this does not apply any transform operation
class MultiLink1DImpl : public INavigationDelegate {
 public:
  std::vector<const DetectorVolume*> dVolumes = {};
  std::vector<ActsScalar> cBoundaries = {};
  BinningValue bValue;

  /// Convenience constructor
  ///
  /// @param detectorVolumes the list of the detector volumes
  /// @param castBoundaries the boundaries inthe cast parameters
  /// @param binningValue the the binning/cast value
  MultiLink1DImpl(const std::vector<const DetectorVolume*>& detectorVolumes,
                  const std::vector<ActsScalar>& castBoundaries,
                  BinningValue binningValue)
      : dVolumes(detectorVolumes),
        cBoundaries(castBoundaries),
        bValue(binningValue) {}

  /// @brief  the function call to be connected with the Delegete
  ///
  /// @param gctx the geometry context for this call
  /// @param position is the 3D global position for this cast
  /// @param direction the direction for for this call
  ///
  /// @note the geometry context is ingored
  /// @note the direction is ignored
  ///
  /// @return the registred volume
  const DetectorVolume* targetVolume(
      [[maybe_unused]] const GeometryContext& ignored, const Vector3& position,
      [[maybe_unused]] const Vector3& direction) const {
    // cast out the value and, find lower bound and return
    ActsScalar cValue = VectorHelpers::cast(position, bValue);
    auto lb = std::lower_bound(cBoundaries.begin(), cBoundaries.end(), cValue);
    size_t b = static_cast<size_t>(std::distance(cBoundaries.begin(), lb) - 1u);

    return dVolumes[b];
  }
};

/// @brief A transformed multi link relation shipt
class TransformedMulitLink1DImpl : public INavigationDelegate {
 public:
  MultiLink1DImpl multiLink;
  Transform3 transform = Transform3::Identity();

  /// Convenience constructor
  ///
  /// @param detectorVolumes the list of the detector volumes
  /// @param castBoundaries the boundaries inthe cast parameters
  /// @param binningValue the the binning/cast value
  /// @param preTransform the transformation before casting
  TransformedMulitLink1DImpl(
      const std::vector<const DetectorVolume*>& detectorVolumes,
      const std::vector<ActsScalar>& castBoundaries, BinningValue binningValue,
      const Transform3& preTransform)
      : multiLink(detectorVolumes, castBoundaries, binningValue),
        transform(preTransform) {}

  /// @brief the function call to be connected with the Delegete
  ///
  /// @param gctx the geometry context for this call
  /// @param position is the 3D global position for this cast to which
  ///                 a transform will be applied first
  /// @param direction the direction for for this call
  ///
  /// @note the geometry context is ingored
  /// @note the direction is ignored
  ///
  /// @return the registred volume
  const DetectorVolume* targetVolume(const GeometryContext& gctx,
                                     const Vector3& position,
                                     const Vector3& direction) const {
    Vector3 tPosition = transform * position;
    return multiLink.targetVolume(gctx, tPosition, direction);
  }
};

/// Helper method to set the outside volume to the internal volumes
///
/// @param volume is the volume to be set
///
/// @note throws an exception if outside link can not be determined
template <typename volume_type>
void setOutsideVolumeLink(volume_type& volume) noexcept(false) {
  // Set the this volume as outside volume of the inserted one
  for (auto& v : volume.volumePtrs()) {
    for (auto& p : v->portalPtrs()) {
      auto singleLinkStored = std::make_shared<SingleLinkImpl>(volume);
      DetectorVolumeLink singleLink;
      singleLink.connect<&SingleLinkImpl::targetVolume>(singleLinkStored.get());
      ManagedDetectorVolumeLink managedLink{std::move(singleLink),
                                            std::move(singleLinkStored)};
      // Get the outside link
      auto [bwd, fwd] = p->volumeLinks();
      if (bwd.implementation == nullptr and fwd.implementation == nullptr) {
        throw std::invalid_argument(
            "DetectorVolumeLinks: outside link can not be determined.");
      }

      NavigationDirection oDir = (bwd.implementation == nullptr)
                                     ? NavigationDirection::Backward
                                     : NavigationDirection::Forward;

      p->updateVolumeLink(oDir, std::move(managedLink));
    }
  }
}

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts

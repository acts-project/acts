// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Experimental/InternalBlueprint.hpp"
#include "Acts/Geometry/Extent.hpp"

#include <exception>
#include <functional>
#include <vector>

namespace Acts {

class DetectorVolume;

using VolumeBuilder = std::function<

/// The volume blueprint
///
/// This class is to be used how volumes are built, they
/// can either be built from a layer structure or from other
/// volumes as a container.
class VolumeBlueprint {
 public:

  /// Constructor for a volume with one layer
  ///
  /// @param gctx The geometry context
  /// @param restriction The extent restriction against which this is checked
  VolumeBlueprint(
    const GeometryContext& gctx,
    const InternalBlueprint& lBlueprint,
    const VolumeBuilder& vBuilder,
    const std::string& name,
    const Extent& restriction = Extent(true)
  )

  /// Constructor for an empty volume
  ///
  /// @param gctx The geometry context
  /// @param restriction The extent restriction against which this is checked
  VolumeBlueprint(
      const GeometryContext& gctx, const Extent& restriction = Extent(true);

};

}  // namespace Acts

/**
namespace Acts {

class DetectorVolume;

/// The void container builder that forwards simply
struct VoidContainerBuilder {
  /// This is the void container builder call, it simply takes
  /// the single volume (required to be one only) from the
  /// provided list and returns it as a the final container.
  ///
  /// @param containerVolumes the volumes from the volume builder
  /// @param name the name to be set to the volume
  ///
  std::shared_ptr<DetectorVolume> operator()(
      std::vector<std::shared_ptr<DetectorVolume>>&& containerVolumes,
      const std::string& name) const noexcept(false) {
    if (containerVolumes.size() != 1 or containerVolumes[0] == nullptr) {
      throw std::invalid_argument(
          "VoidContainerBuilder: exacly one Volume has to be provided.");
    }
    // The single volume is the new container
    auto volume = containerVolumes[0];
    volume->setName(name);
    return volume;
  }
};

/// This defines how the sub volumes are build
///
///
using VolumeBuilder =
    std::function<std::vector<std::shared_ptr<DetectorVolume>>(
        const GeometryContext&,
        const Extent& extent,
        const std::vector<InternalBlueprint>& InternalBlueprints,
        std::string&)>;

/// This defines how the container is build
///
/// In case the volume builder returns only one volume, this is the
/// container volume already;
using ContainerBuilder = std::function<std::shared_ptr<DetectorVolume>(
    std::vector<std::shared_ptr<DetectorVolume>>&& containerVolumes,
    const std::string&)>;

/// The volume blue print that is used to describe the building
/// instruction of a volume that contains layer volumes
///
class VolumeBlueprint {
 public:
  /// Constructor from arguments for a volume
  ///
  /// @param extent is the volume extent
  /// @param InternalBlueprints is the list of layer blueprints
  /// @param vBuilder is the function for layer volume building
  /// @param cBuilder is the function for container volume building
  VolumeBlueprint(const Extent& extent,
                  const std::vector<InternalBlueprint>& InternalBlueprints,
                  VolumeBuilder vBuilder,
                  ContainerBuilder cBuilder = VoidContainerBuilder());

  VolumeBlueprint() = delete;

  /// @return the current extent - const access
  const Extent& extent() const;

  /// access to the current extent - non-const access
  Extent& extent();

  /// @return the layer blue prints
  const std::vector<InternalBlueprint>& InternalBlueprints() const;

  /// @return the volume builder function
  const VolumeBuilder& volumeBuilder() const;

  /// @return the container builder function
  const ContainerBuilder& containerBuilder() const;

 private:
  /// The extent for this volume
  Extent m_extent;

  /// The contained layer blueprints
  std::vector<InternalBlueprint> m_InternalBlueprints;

  /// The volume builder to be used with this blue print
  VolumeBuilder m_volumeBuilder;

  /// Teh container builder used with the blue print
  ContainerBuilder m_containerBuilder;
};

inline const Extent& VolumeBlueprint::extent() const {
  return m_extent;
}

inline Extent& VolumeBlueprint::extent() {
  return m_extent;
}

inline const std::vector<InternalBlueprint>& VolumeBlueprint::InternalBlueprints()
    const {
  return m_InternalBlueprints;
}

inline const VolumeBuilder& VolumeBlueprint::volumeBuilder()
    const {
  return m_volumeBuilder;
}

inline const ContainerBuilder&
VolumeBlueprint::containerBuilder() const {
  return m_containerBuilder;
}

}  // namespace Acts
*/

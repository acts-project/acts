// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Detector/interface/IDetectorComponentBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <memory>
#include <string>

namespace Acts {
class ISurfaceMaterial;
}

namespace ActsTests {

/// @brief A mockup volume builder, it generates volumes with
/// a single surface filled in in order to use the CylindricalContainerBuilder
/// infrastructure.
template <typename surface_type, typename surface_bounds_type>
class CylindricalVolumeBuilder
    : public Acts::Experimental::IDetectorComponentBuilder {
 public:
  /// @brief Constructor from arguments
  /// @param transform the positioning of the volume
  /// @param vBounds the volume bounds
  /// @param sBounds the surface bounds
  /// @param vName the volume name
  /// @param material the surface material
  CylindricalVolumeBuilder(
      const Acts::Transform3& transform,
      const Acts::CylinderVolumeBounds& vBounds,
      const surface_bounds_type& sBounds, const std::string& vName,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr)
      : IDetectorComponentBuilder(),
        m_transform(transform),
        m_volumeBounds(vBounds),
        m_surfaceBounds(sBounds),
        m_name(vName),
        m_material(std::move(material)) {}

  /// @brief Construct the detector component
  /// @param gctx the geometry context at construction time
  Acts::Experimental::DetectorComponent construct(
      [[maybe_unused]] const Acts::GeometryContext& gctx) const final {
    // The outgoing root volumes
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>
        rootVolumes;

    // Ingredients
    auto surface = Acts::Surface::makeShared<surface_type>(
        (m_transform), std::make_shared<surface_bounds_type>(m_surfaceBounds));
    surface->assignSurfaceMaterial(m_material);

    auto bounds = std::make_unique<Acts::CylinderVolumeBounds>(m_volumeBounds);
    auto portalGenerator = Acts::Experimental::defaultPortalGenerator();
    auto volume = Acts::Experimental::DetectorVolumeFactory::construct(
        portalGenerator, gctx, m_name, m_transform, std::move(bounds),
        {surface}, {}, Acts::Experimental::tryNoVolumes(),
        Acts::Experimental::tryAllPortalsAndSurfaces());

    // Add to the roots
    rootVolumes.push_back(volume);

    Acts::Experimental::DetectorComponent::PortalContainer dContainer;
    for (auto [ip, p] : Acts::enumerate(volume->portalPtrs())) {
      dContainer[ip] = p;
    }
    return Acts::Experimental::DetectorComponent{
        {volume},
        dContainer,
        Acts::Experimental::RootDetectorVolumes{
            rootVolumes, Acts::Experimental::tryRootVolumes()}};
  }

 private:
  Acts::Transform3 m_transform;
  Acts::CylinderVolumeBounds m_volumeBounds;
  surface_bounds_type m_surfaceBounds;
  std::string m_name;
  std::shared_ptr<const Acts::ISurfaceMaterial> m_material = nullptr;
};

/// @brief A mockup container builder, it generates a container with
/// serval cylindrical volumes in it.
std::shared_ptr<const Acts::Experimental::Detector> buildCylindricalDetector(
    const Acts::GeometryContext& tContext);

}  // namespace ActsTests

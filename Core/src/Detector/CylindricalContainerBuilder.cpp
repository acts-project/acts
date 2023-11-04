// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/CylindricalContainerBuilder.hpp"

#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/DetectorVolumeBuilder.hpp"
#include "Acts/Detector/VolumeStructureBuilder.hpp"
#include "Acts/Detector/detail/CylindricalDetectorHelper.hpp"
#include "Acts/Detector/interface/IGeometryIdGenerator.hpp"
#include "Acts/Detector/interface/IRootVolumeFinderBuilder.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"

#include <algorithm>
#include <ostream>
#include <stdexcept>
#include <utility>

namespace Acts {
namespace Experimental {
class DetectorVolume;
}  // namespace Experimental
}  // namespace Acts

namespace {

/// @brief Helper method to connect/wrap volumes or containers
///
/// @tparam object_collection either a vector of volumes or containers
/// @param gctx The geometry context of the this call
/// @param objects The object vector
/// @param binnning the chosen binning
/// @param logLevel the logging output level
///
/// @note no checking on consistency is done as the caller container builder
/// is already checked at construction
///
/// @return a newly built container
template <typename object_collection>
Acts::Experimental::DetectorComponent::PortalContainer connect(
    const Acts::GeometryContext& gctx, object_collection& objects,
    const std::vector<Acts::BinningValue>& binning,
    Acts::Logging::Level logLevel) {
  // Return container object
  Acts::Experimental::DetectorComponent::PortalContainer rContainer;
  if (binning.size() == 1u) {
    Acts::BinningValue bv = binning.front();
    // 1-dimensional binning options
    switch (bv) {
      case Acts::binR: {
        rContainer =
            Acts::Experimental::detail::CylindricalDetectorHelper::connectInR(
                gctx, objects, {}, logLevel);
      } break;
      case Acts::binZ: {
        rContainer =
            Acts::Experimental::detail::CylindricalDetectorHelper::connectInZ(
                gctx, objects, {}, logLevel);
      } break;
      case Acts::binPhi: {
        rContainer =
            Acts::Experimental::detail::CylindricalDetectorHelper::connectInPhi(
                gctx, objects, {}, logLevel);
      } break;
      default:
        break;
    }
  } else if (binning ==
                 std::vector<Acts::BinningValue>{Acts::binZ, Acts::binR} and
             objects.size() == 2u) {
    rContainer =
        Acts::Experimental::detail::CylindricalDetectorHelper::wrapInZR(
            gctx, objects, logLevel);
  }
  return rContainer;
}
}  // namespace

Acts::Experimental::CylindricalContainerBuilder::CylindricalContainerBuilder(
    const Acts::Experimental::CylindricalContainerBuilder::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : IDetectorComponentBuilder(), m_cfg(cfg), m_logger(std::move(logger)) {
  // Check if builders are present
  if (m_cfg.builders.empty()) {
    throw std::invalid_argument(
        "CylindricalContainerBuilder: no sub builders provided.");
  }
  // Check if binning value is correctly chosen
  if (m_cfg.binning.size() == 1u) {
    // 1-dimensional case
    auto b = m_cfg.binning.front();
    if (b != Acts::binR and b != Acts::binZ and b != Acts::binPhi) {
      throw std::invalid_argument(
          "CylindricalContainerBuilder: 1D binning only supported in z, r, or "
          "phi");
    }
  } else if (m_cfg.binning.size() == 2u) {
    // 2-dimensional case, this is for wrapping
    if (m_cfg.binning !=
        std::vector<Acts::BinningValue>{Acts::binZ, Acts::binR}) {
      throw std::invalid_argument(
          "CylindricalContainerBuilder: 2D binning only supports wrapping in "
          "z-r.");
    } else if (m_cfg.builders.size() != 2u) {
      // Wrapping needs exactly one inner (volume or container) and one outer
      // volume
      throw std::invalid_argument(
          "CylindricalContainerBuilder: 2D wrapping in z-r requires exactly "
          "two builders.");
    }
  }
}

Acts::Experimental::CylindricalContainerBuilder::CylindricalContainerBuilder(
    const Acts::Experimental::Blueprint::Node& bpNode,
    Acts::Logging::Level logLevel)
    : IDetectorComponentBuilder(),
      m_logger(getDefaultLogger(bpNode.name + "_cont", logLevel)) {
  if (bpNode.boundsType != VolumeBounds::BoundsType::eCylinder) {
    throw std::invalid_argument(
        "CylindricalContainerBuilder: boundary type must be cylinder - for "
        "building from a blueprint node.");
  }

  std::vector<std::shared_ptr<const IDetectorComponentBuilder>> builders;
  for (const auto& child : bpNode.children) {
    if (child->isLeaf()) {
      // Volume structure
      VolumeStructureBuilder::Config vsCfg;
      vsCfg.transform = child->transform;
      vsCfg.boundsType = child->boundsType;
      vsCfg.boundValues = child->boundaryValues;
      vsCfg.auxiliary = "*** acts auto-generated shape builder ***";
      auto vsBuilder = std::make_shared<VolumeStructureBuilder>(
          vsCfg, getDefaultLogger(child->name + "_shape", logLevel));
      // Detector volume builder
      DetectorVolumeBuilder::Config dvCfg;
      dvCfg.name = child->name;
      dvCfg.externalsBuilder = vsBuilder;
      dvCfg.internalsBuilder = child->internalsBuilder;
      dvCfg.auxiliary = "*** acts auto-generated volume builder ***";
      // Add the builder
      m_cfg.builders.push_back(std::make_shared<DetectorVolumeBuilder>(
          dvCfg, getDefaultLogger(child->name, logLevel)));
    } else {
      // This evokes the recursive stepping down the tree
      m_cfg.builders.push_back(
          std::make_shared<CylindricalContainerBuilder>(*child, logLevel));
    }
  }

  m_cfg.binning = bpNode.binning;
  m_cfg.auxiliary = "*** acts auto-generated from proxy ***";
  m_cfg.geoIdGenerator = bpNode.geoIdGenerator;
  m_cfg.rootVolumeFinderBuilder = bpNode.rootVolumeFinderBuilder;
}

Acts::Experimental::DetectorComponent
Acts::Experimental::CylindricalContainerBuilder::construct(
    const GeometryContext& gctx) const {
  // Return container object
  DetectorComponent::PortalContainer rContainer;
  bool atNavigationLevel = true;

  // Create the indivudal components, collect for both outcomes
  std::vector<DetectorComponent> components;
  ACTS_DEBUG("Building container from " << m_cfg.builders.size()
                                        << " components.");
  // Check through the component volumes - if every builder only
  // built exactly one volume, you are at pure navigation level
  // Collect the volumes
  std::vector<std::shared_ptr<DetectorVolume>> volumes;
  std::vector<DetectorComponent::PortalContainer> containers;
  std::vector<std::shared_ptr<DetectorVolume>> rootVolumes;
  // Run through the builders
  std::for_each(
      m_cfg.builders.begin(), m_cfg.builders.end(), [&](const auto& builder) {
        auto [cVolumes, cContainer, cRoots] = builder->construct(gctx);
        atNavigationLevel = (atNavigationLevel and cVolumes.size() == 1u);
        // Collect individual components, volumes, containers, roots
        volumes.insert(volumes.end(), cVolumes.begin(), cVolumes.end());
        containers.push_back(cContainer);
        rootVolumes.insert(rootVolumes.end(), cRoots.volumes.begin(),
                           cRoots.volumes.end());
      });
  // Navigation level detected, connect volumes (cleaner and faster than
  // connect containers)
  if (atNavigationLevel) {
    ACTS_VERBOSE(
        "Component volumes are at navigation level: connecting volumes.");
    // Connect volumes
    rContainer = connect(gctx, volumes, m_cfg.binning, logger().level());
  } else {
    ACTS_VERBOSE("Components contain sub containers: connect containers.");
    // Connect containers
    rContainer = connect(gctx, containers, m_cfg.binning, logger().level());
  }
  ACTS_VERBOSE("Number of root volumes: " << rootVolumes.size());

  // Check if a root volume finder is provided
  if (m_cfg.rootVolumeFinderBuilder) {
    // Return the container
    return Acts::Experimental::DetectorComponent{
        {},
        rContainer,
        RootDetectorVolumes{
            rootVolumes,
            m_cfg.rootVolumeFinderBuilder->construct(gctx, rootVolumes)}};
  }

  // Geometry Id generation
  if (m_cfg.geoIdGenerator != nullptr) {
    ACTS_DEBUG("Assigning geometry ids to the detector");
    auto cache = m_cfg.geoIdGenerator->generateCache();
    if (m_cfg.geoIdReverseGen) {
      std::for_each(rootVolumes.rbegin(), rootVolumes.rend(), [&](auto& v) {
        m_cfg.geoIdGenerator->assignGeometryId(cache, *v);
      });
    } else {
      std::for_each(rootVolumes.begin(), rootVolumes.end(), [&](auto& v) {
        m_cfg.geoIdGenerator->assignGeometryId(cache, *v);
      });
    }
  }

  // Return the container
  return Acts::Experimental::DetectorComponent{
      {}, rContainer, RootDetectorVolumes{rootVolumes, tryRootVolumes()}};
}

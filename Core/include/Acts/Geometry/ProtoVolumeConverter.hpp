// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/DetectorHelper.hpp"
#include "Acts/Geometry/DetectorVolume.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationDelegates.hpp"
#include "Acts/Geometry/ProtoDetector.hpp"
#include "Acts/Geometry/detail/NavigationStateUpdators.hpp"
#include "Acts/Geometry/detail/PortalGenerators.hpp"
#include "Acts/Geometry/detail/SurfaceCandidatesUpdators.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <tuple>

namespace Acts {
namespace Experimental {

class Portal;

using DetectorVolumeExternals =
    std::tuple<Transform3, std::unique_ptr<CylinderVolumeBounds>>;

// Struct to convert ProtoVolumes into a cylincrical DetectorVolumes
struct ConcentricCylinderConverter {
  /// The proto volume that needs to be converted
  ProtoVolume protoVolume;

  ///  Create cylindrical volume bounds
  ///
  /// @param gctx The geometry context
  DetectorVolumeExternals create(const GeometryContext& /*ignored*/) {
    // Get the extent of this volume and translate
    const auto& pvExtent = protoVolume.extent;
    ActsScalar z = pvExtent.medium(binZ);
    Transform3 transform = Transform3::Identity();
    transform.pretranslate(Vector3(0., 0., z));
    // Now the shape
    ActsScalar rI = pvExtent.min(binR);
    ActsScalar rO = pvExtent.max(binR);
    ActsScalar hZ = 0.5 * pvExtent.span(binZ);
    ActsScalar hPhi = M_PI;
    ActsScalar aPhi = 0.;
    if (pvExtent.constrains(binPhi)) {
      hPhi = 0.5 * pvExtent.span(binPhi);
      aPhi = pvExtent.medium(binPhi);
    }
    // Return a tuple of transform and bounds
    auto bounds =
        std::make_unique<CylinderVolumeBounds>(rI, rO, hZ, hPhi, aPhi);
    DetectorVolumeExternals rTuple = {std::move(transform), std::move(bounds)};
    return rTuple;
  }
};

using DetectorVolumeInternals =
    std::tuple<std::vector<std::shared_ptr<Surface>>,
               std::vector<std::shared_ptr<DetectorVolume>>,
               SurfaceCandidatesUpdator>;

/// A struct that creates empty internals of a volume
struct EmptyInternals {
  /// The proto volume that needs to be converted, ignored here
  ProtoVolume protoVolume;

  /// @brief Create an internal structure for only portals
  ///
  /// @return a tuple of surfaces, volumes, updators
  DetectorVolumeInternals create(const GeometryContext& /* ingored */) {
    std::vector<std::shared_ptr<Surface>> noSurfaces = {};
    std::vector<std::shared_ptr<DetectorVolume>> noVolumes = {};
    auto allPortals = detail::allPortals();
    DetectorVolumeInternals rTuple = {
        std::move(noSurfaces), std::move(noVolumes), std::move(allPortals)};
    return rTuple;
  }
};

/// A struct that makes a default portal generator
struct DefaultPortalsConverter {
  /// The proto volume that needs to be converted, ignored here
  ProtoVolume protoVolume;

  /// This method creates a standard portal generator
  auto create(const GeometryContext& /* ingored */) {
    return detail::defaultPortalGenerator();
  }
};

/// Definition of a shell builder, it builds a detector volume and
/// applies its shell to it
template <typename Volume = ConcentricCylinderConverter,
          typename Portals = DefaultPortalsConverter,
          typename Internals = EmptyInternals>
struct SingleBlockBuilder {
  /// The proto volume that needs to be converted
  ProtoVolume protoVolume;

  /// Convert a proto volume into a detector volume
  ///
  /// @tparam Volume how the volume is handled (bounds, position)
  /// @tparam Portals how the portals are handled
  /// @tparam InternalsHandling how the internals are handled
  ///
  /// @param shell The input shell
  /// @param gctx The geometry context
  /// @param logLevel is a screen output log level
  ///
  /// @return a newly created DetectorVolume
  void operator()(DetectorBlock& dBlock, const GeometryContext& gctx,
                  Acts::Logging::Level logLevel = Acts::Logging::INFO) {
    // Screen output
    ACTS_LOCAL_LOGGER(getDefaultLogger(
        "SingleBlockBuilder   [ " + protoVolume.name + " ]", logLevel));

    ACTS_DEBUG("Building single volume '" << protoVolume.name << "'.");
    auto& dVolumes = std::get<DetectorVolumes>(dBlock);
    auto& dContainer = std::get<ProtoContainer>(dBlock);
    // Externals
    auto [transform, bounds] = Volume{protoVolume}.create(gctx);
    // Internals
    auto [surfaces, volumes, updator] = Internals{protoVolume}.create(gctx);
    // Portals
    auto portals = Portals{protoVolume}.create(gctx);
    // Construct the detector volume
    auto dVolume = DetectorVolumeFactory::construct(
        portals, gctx, protoVolume.name, transform, std::move(bounds), surfaces,
        volumes, std::move(updator));
    // Update the detector block
    dVolumes.push_back(dVolume);
    for (auto [ip, p] : enumerate(dVolume->portalPtrs())) {
      ACTS_VERBOSE(" - adding portal " << ip << " to the proto container.");
      dContainer[ip] = p;
    }
    ACTS_VERBOSE(" - total number of portals added: " << dContainer.size());
  }
};

struct ContainerBlockBuilder {
  /// The proto volume that needs to be converted
  ProtoVolume protoVolume;

  /// Convert a proto volume into a detector volume
  ///
  /// @tparam Volume how the volume is handled (bounds, position)
  /// @tparam Portals how the portals are handled
  /// @tparam InternalsHandling how the internals are handled
  ///
  /// @param dBlock The detector block to be built
  /// @param gctx The geometry context
  /// @param logLevel is a sreen output log level
  ///
  /// @return a newly created DetectorVolume
  void operator()(DetectorBlock& dBlock, const GeometryContext& gctx,
                  Acts::Logging::Level logLevel = Acts::Logging::INFO) {
    // Screen output
    ACTS_LOCAL_LOGGER(getDefaultLogger(
        "ContainerBlockBuilder [ " + protoVolume.name + " ]", logLevel));

    ACTS_DEBUG("Building container volume '" << protoVolume.name << "'.");
    auto& dVolumes = std::get<DetectorVolumes>(dBlock);
    // Only one dimensional connection so far possible
    if (protoVolume.constituentBinning.size() == 1u) {
      // Get the binning value
      BinningValue bValue = protoVolume.constituentBinning.front().binvalue;
      // Run the block builders
      std::vector<ProtoContainer> dContainers = {};
      ACTS_VERBOSE(" - this container has "
                  << protoVolume.constituentVolumes.size() << " constituents");
      for (auto& cv : protoVolume.constituentVolumes) {
        // Collect the constituent volumes
        DetectorBlock cBlock;
        cv.blockBuilder(cBlock, gctx, logLevel);
        // Collect the detector volumes & container of the constituent
        const auto& cVolumes = std::get<DetectorVolumes>(cBlock);
        auto& cContainer = std::get<ProtoContainer>(cBlock);
        dVolumes.insert(dVolumes.end(), cVolumes.begin(), cVolumes.end());
        // Collect the proto containers
        dContainers.push_back(cContainer);
      }
      auto& dContainer = std::get<ProtoContainer>(dBlock);
      dContainer = connectContainers(gctx, bValue, dContainers, {}, logLevel);
    }

    // throw exception
    throw std::invalid_argument(
        "ContainerBlockBuilder: no binning value provided.");
  }
};

}  // namespace Experimental

}  // namespace Acts

// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/ProtoBinning.hpp"
#include "Acts/Detector/ProtoSupport.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"
#include "Acts/Detector/interface/ISurfacesProvider.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace Acts {
namespace Experimental {

/// @brief This is a builder of layer structures to be contained
/// within a DetectorVolume, it extends the IInternalStructureBuilder
/// interface and provides the internal structure components of
/// DetectorVolume objects to be constructed.
///
/// It uses the IndexedSurfaceGrid to bin the internal surfaces,
/// and allows for additional support surfaces that are added to the
/// structure and indexing mechanism. Those support structures can
/// also be approximated by planar surfaces, in order to facilitate
/// vectorization of surface intersection calls.
///
/// The binning can be chosen with a so called `expansion`, a number
/// which indicates the configured expanded bin window in which the
/// surfaces are going to be filled, the details to this strategy
/// can be found in the IndexedGridFiller and IndexedSurfacesGenerator
/// classes.
///
/// No sub volumes are added to this structure builders, hence,
/// the DetectorVolumeFinder navigation delegate uses the "NoopFinder"
/// breakpoint to indicate the bottom of the volume hierarchy.
///
class LayerStructureBuilder : public IInternalStructureBuilder {
 public:
  /// @brief A holder struct for surfaces
  class SurfacesHolder final : public ISurfacesProvider {
   public:
    /// Constructor with predefined surfaces
    /// @param isurfaces is the vector of surfaces
    SurfacesHolder(std::vector<std::shared_ptr<Surface>> isurfaces)
        : m_surfaces(std::move(isurfaces)) {}

    /// Return the surfaces from the holder
    /// @param gctx is the geometry context
    std::vector<std::shared_ptr<Surface>> surfaces(
        [[maybe_unused]] const GeometryContext& gctx) const final {
      return m_surfaces;
    }

   private:
    std::vector<std::shared_ptr<Surface>> m_surfaces = {};
  };

  /// @brief Configuration struct for the LayerStructureBuilder
  ///
  /// It contain:
  /// - a source of the surfaces to be built
  /// - a definition of surface binning on this layer
  /// - a definition of supports to be built
  struct Config {
    /// Connection point for a function to provide surfaces
    std::shared_ptr<ISurfacesProvider> surfacesProvider = nullptr;
    /// Definition of Supports
    std::vector<ProtoSupport> supports = {};
    /// Definition of Binnings
    std::vector<ProtoBinning> binnings = {};
    /// Polyhedron approximations
    unsigned int nSegments = 1u;
    /// Extra information, mainly for screen output
    std::string auxiliary = "";
  };

  /// Constructor
  ///
  /// @param cfg is the configuration struct
  /// @param logger logging instance for screen output
  LayerStructureBuilder(const Config& cfg,
                        std::unique_ptr<const Logger> logger = getDefaultLogger(
                            "LayerStructureBuilder", Logging::INFO));

  /// The interface definition for internal structure creation
  ///
  /// @param gctx the geometry context at the creation of the internal structure
  ///
  /// @return a consistent set of detector volume internals
  InternalStructure construct(const GeometryContext& gctx) const final;

 private:
  /// configuration object
  Config m_cfg;

  /// Private access method to the logger
  const Logger& logger() const { return *m_logger; }

  /// logging instance
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Experimental
}  // namespace Acts

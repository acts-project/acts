// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <functional>
#include <optional>

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
  /// @brief Support parameter defintions
  struct Support {
    /// Define whether you want to build support structures
    std::array<ActsScalar, 5u> values = {};
    /// The surface type to be built
    Surface::SurfaceType type = Surface::SurfaceType::Other;
    /// Define in which values the support should be constrained
    std::vector<BinningValue> constraints = s_binningValues;
    /// Potential splits into planar approximations
    unsigned int splits = 1u;
    /// The (optional) layer transform
    std::optional<Transform3> transform = std::nullopt;
  };

  /// @brief The surface binning definition
  struct Binning {
    /// Define the binning of the surfaces
    BinningData data;
    /// An expansion for the filling (in bins)
    size_t expansion = 0u;
  };

  /// @brief Configuration struct for the LayerStructureBuilder
  ///
  /// It contain:
  /// - a source of the surfaces to be built
  /// - a definition of surface binning on this layer
  /// - a definition of supports to be built
  struct Config {
    /// Connection point for a function to provide surfaces
    std::function<std::vector<std::shared_ptr<Surface>>()> surfaces;
    /// Definition of Supports
    std::vector<Support> supports = {};
    /// Definition of Binnings
    std::vector<Binning> binnings = {};
    /// Polyhedron approximations
    unsigned int nSegments = 1u;
    /// Extra information, mainly for screen output
    std::string auxilliary = "";
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

  /// Private acces method to the logger
  const Logger& logger() const { return *m_logger; }

  /// logging instance
  std::unique_ptr<const Logger> m_logger;
};

}  // namespace Experimental
}  // namespace Acts

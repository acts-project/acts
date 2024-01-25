// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/interface/ISurfacesProvider.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorSurfaceFactory.hpp"
#include "Acts/Detector/KdtSurfacesProvider.hpp" 

#include "G4GDMLParser.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

namespace Acts {
namespace Experimental {

/// @brief A surface provider that extracts surfaces from a gdml file
/// 
/// This provider extracts volumes from a gdml file based on 
/// the names and returns them as a vector of shared pointers 
/// to surfaces. 
///
/// Additionally, it can be configured to return a range-based 
/// subset of all the surfaces selected by the name.
/// The latter is achieved by using a KdtSurfacesProvider 
/// internally.
/// 
/// @note If the names of the volumes to be extracted are unique,
/// and the range-based selection is not needed, instantialization of
/// the KdtSurfacesProvider can be skipped. In this case, the 
/// template parameters do not affect the behaviour of the provider.
///
/// @tparam kDim The number of dimensions for the KDTree
/// @tparam bSize The maximum number of surfaces per KDTree leaf
/// @tparam reference_generator The reference generator for the KDTree
template <std::size_t kDim = 2u, std::size_t bSize = 100u,
  typename reference_generator =
    detail::PolyhedronReferenceGenerator<1u, false>>
class Geant4SurfaceProvider : public Acts::Experimental::ISurfacesProvider {
 public:
  /// Nested configuration struct
  struct Config {
    /// The path of the gdml file
    std::string gdmlPath = "";

    /// The names of the surfaces to
    /// extract from the gdml file
    std::vector<std::string> surfaceNames;

    /// Whether to match the names exactly
    bool exactMatch = false;
  };

  /// Optional configuration for the KdtSurfaces
  struct kdtOptions {
    /// A set of ranges to separate the surfaces
    Acts::RangeXD<kDim,Acts::ActsScalar> range;

    std::array<Acts::BinningValue,kDim> binningValues;

    std::size_t leafSize = bSize;
    
    reference_generator rgen;
  };

  /// Constructor
  ///@param config The configuration struct
  Geant4SurfaceProvider(const Config& config,
    const kdtOptions& options = kdtOptions()) {
      if (config.gdmlPath.empty()) {
        throw std::invalid_argument("Geant4SurfaceProvider: no gdml file provided");
      }
      if (config.surfaceNames.empty()) {
        throw std::invalid_argument(
            "Geant4SurfaceProvider: no surface names provided");
      }
  
      m_cfg = config;
      m_kdtOptions = options;

      /// Read the gdml file and get the world volume
      G4GDMLParser parser;
      parser.Read(m_cfg.gdmlPath);
      m_g4World = parser.GetWorldVolume();
  
      if (m_g4World == nullptr) {
        throw std::invalid_argument(
            "Geant4SurfaceProvider: No g4World initialized");
      }
    };

  /// Destructor
  ~Geant4SurfaceProvider() = default;

  std::vector<std::shared_ptr<Acts::Surface>> surfaces(
    [[maybe_unused]] const Acts::GeometryContext& gctx) const override {
      /// Surface selector to look in the G4 world
      auto surfaceSelector =
          std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(
              m_cfg.surfaceNames, m_cfg.exactMatch);

      /// Surface factory options
      Acts::Geant4DetectorSurfaceFactory::Options g4SurfaceOptions;

      /// It doesn't matter which selector is used here
      g4SurfaceOptions.passiveSurfaceSelector = surfaceSelector;
    
      /// Generate the surface cache
      Acts::Geant4DetectorSurfaceFactory::Cache g4SurfaceCache;
      G4Transform3D g4ToWorld;

      /// Find and store surfaces in the cache object
      Acts::Geant4DetectorSurfaceFactory{}.construct(
        g4SurfaceCache, g4ToWorld, *m_g4World, g4SurfaceOptions);

      auto surfaces = g4SurfaceCache.passiveSurfaces;

      /// If range is not set, return all surfaces
      if (m_kdtOptions.range.degenerate()) {
        return surfaces;
      }

      /// Otherwise, select the surfaces based on the range
      auto kdtSurfaces = 
        Acts::Experimental::KdtSurfaces<kDim,bSize,reference_generator>
          (gctx,surfaces,m_kdtOptions.binningValues);

      return kdtSurfaces.surfaces(m_kdtOptions.range);
  };

 private:
  Config m_cfg;

  kdtOptions m_kdtOptions;

  G4VPhysicalVolume* m_g4World = nullptr;
};

}  // namespace Experimental
}  // namespace Acts

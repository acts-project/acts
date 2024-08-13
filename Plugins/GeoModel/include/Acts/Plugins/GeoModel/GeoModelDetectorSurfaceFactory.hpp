// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Plugins/GeoModel/IGeoShapeConverter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <GeoModelHelpers/getChildNodesWithTrf.h>

class GeoShape;
class GeoFullPhysVol;

namespace Acts {

struct GeoModelTree;
class Surface;

/// A factory to convert GeoModel volume into sensitive
/// or passive surfaces which are filled into a Cache object,
/// also the create GeoModelDetectorElements which are also
/// returned.
class GeoModelDetectorSurfaceFactory {
 public:
  /// Collect the passive surfaces, bool whether it should be
  /// added as an "always try, i.e. assignToAll=true" surface
  using GeoModelPassiveSurface = std::tuple<std::shared_ptr<Surface>, bool>;

  // Configuration struct for the Detector surface factory
  struct Config {
    // /// List for names to match
    std::vector<std::string> nameList;
    /// List for materials to match
    std::vector<std::string> materialList;

    /// boolean flag to build subvolumes
    bool convertSubVolumes = false;
  };

  /// Nested cache that records the conversion status
  struct Cache {
    /// The created detector elements and their surfaces
    std::vector<GeoModelSensitiveSurface> sensitiveSurfaces;
    /// The created passive representation surfaces
    std::vector<GeoModelPassiveSurface> passiveSurfaces;
  };

  /// The options to steer the conversion
  struct Options {
    std::vector<std::string> queries = {};
  };

  /// The GeoModel detector element factory
  ///
  /// @param cfg the configuration struct
  /// @param mlogger a screen output logger
  GeoModelDetectorSurfaceFactory(
      const Config& cfg,
      std::unique_ptr<const Logger> mlogger = getDefaultLogger(
          "GeoModelDetectorSurfaceFactory", Acts::Logging::WARNING));

  /// Construction method of the detector elements
  ///
  /// @param cache [in,out] into which the Elements are filled
  /// @param gctx the geometry context
  /// @param geoModelTree the gjeo model tree
  /// @param options to steer the conversion
  ///
  /// @note this method will call the recursive construction
  void construct(Cache& cache, const GeometryContext& gctx,
                 const GeoModelTree& geoModelTree, const Options& options);

 private:
  /// The configuration struct
  Config m_cfg;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to the logger
  const Logger& logger() const { return *m_logger; }

  /// Private helper method that does the conversion to surfaces
  ///
  /// @param geoPV the pointer to the physical volume
  /// @param transform the transform to the global frame
  /// @param sensitives the vector which we fill with the converted surfaces
  void convertSensitive(
      PVConstLink geoPV, const Acts::Transform3& transform,
      std::vector<Acts::GeoModelSensitiveSurface>& sensitives);

  /// Private helper method that finds the subvolumes from the top-level full
  /// physical volume
  ///
  /// @param geoPV the pointer to the physical volume
  /// @param matchFunc the function that checks if the volume matches the name and material requirements
  std::vector<GeoChildNodeWithTrf> findAllSubVolumes(
      PVConstLink geoPV,
      std::function<bool(std::string, PVConstLink)> matchFunc);
};

}  // namespace Acts

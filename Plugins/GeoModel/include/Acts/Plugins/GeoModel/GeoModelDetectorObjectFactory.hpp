// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Plugins/GeoModel/GeoModelToDetectorVolume.hpp"
#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/BoundFactory.hpp"

#include <GeoModelHelpers/getChildNodesWithTrf.h>

#include "GeoModelKernel/GeoDefinitions.h"

class GeoShape;
struct GeoModelTree;
class Surface;
namespace Acts {
/** @brief Factory class to convert GeoModel objects into Acts volumes. Currently,  */
class GeoModelDetectorObjectFactory {
 public:
  using GeoModelBoundingBox = std::shared_ptr<Experimental::DetectorVolume>;

  struct Options {
    std::vector<std::string> queries = {};
  };

  // substring matching for additional parameters
  // empty vectors will be matched to everything
  struct Config {
    // List for names to match
    std::vector<std::string> nameList;
    // List for materials to match
    std::vector<std::string> materialList;

    // boolean flag to build subvolumes
    bool convertSubVolumes = false;

    // flag to build the desired bounding Boxes
    std::vector<std::string> convertBox;
  };

  struct Cache {
    // The created detector elements and their surfaces
    std::vector<GeoModelSensitiveSurface> sensitiveSurfaces;
    // The created representation of bounding box
    std::vector<GeoModelBoundingBox> boundingBoxes;
    /** @brief Pointer to the surface bound factory  */
    std::shared_ptr<SurfaceBoundFactory> surfBoundFactory = std::make_shared<SurfaceBoundFactory>();
    /** @brief Pointer to the volume bound factory */
    std::shared_ptr<VolumeBoundFactory> volumeBoundFactory = std::make_shared<VolumeBoundFactory>();
  };

  explicit GeoModelDetectorObjectFactory(
      const Config& cfg,
      std::unique_ptr<const Logger> mlogger = getDefaultLogger(
          "GeoModelDetectorObjectFactory", Acts::Logging::WARNING));

  void construct(Cache& cache, const GeometryContext& gctx,
                 const GeoModelTree& geoModelTree, const Options& options);


 private:

  /** @brief Convert the GeoPhysVol into a sensitive Acts::Surface.
   *  @param geoPV: Pointer to the GeoPhysVol to convert
   *  @param transform: Placement of the resulting surface in the world
   *  @param boundFactory: Reference to the BoundFactory to avoid duplicated bounds 
   *                       across similar surfaces
   * @param sensitives: Output vector into which the new converted surface is pushed */
  void convertSensitive(const PVConstLink& geoPV,
                        const Acts::Transform3& transform,
                        SurfaceBoundFactory& boundFactory,
                        std::vector<GeoModelSensitiveSurface>& sensitives);
  /** @brief Find all sub volumes of a passed volume that are 
   *         good for sensitive detector conversion
   *  @param vol: Pointer to the GeoPhysVol to search through
   *  @return A vector of GeoChildNodeWithTrf containing the information about the
   *          volumes to convert and their placement w.r.t. the passed volume */
  std::vector<GeoChildNodeWithTrf> findAllSubVolumes(const PVConstLink& vol) const;
  /** @brief Checks whether the volume name satisfies the user-defined tokens and/or
   *         the material of physical volume does it.
   *  @param name: Name of the physical volume to test. Usually, it's the GeoNameTag or 
   *               the published GeoFullPhysVol entry
   *  @param physvol: Reference to the PhysicalVolume to additionally check material compatibility */
  bool matches(const std::string& name, const PVConstLink& physvol) const;
  
  void convertFpv(const std::string& name, const GeoFullPhysVol* fpv, Cache& cache,
    const GeometryContext& gctx);
bool convertBox(const std::string& name);



  std::unique_ptr<const Logger> m_logger;
  Config m_cfg;

  const Logger& logger() const { return *m_logger; }
};
}  // namespace Acts

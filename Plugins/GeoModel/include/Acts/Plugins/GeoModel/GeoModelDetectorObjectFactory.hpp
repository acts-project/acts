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
#include "Acts/Utilities/BoundFactory.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <GeoModelHelpers/getChildNodesWithTrf.h>

#include "GeoModelKernel/GeoDefinitions.h"

class GeoShape;
struct GeoModelTree;
class Surface;

namespace Acts {

/// @brief Factory class to convert GeoModel objects into Acts volumes and surfaces. The surface conversion
///        process constructs surfaces from GeoTrd, GeoBox, GeoTube,
///        GeoSimplePolygonBrep volume shapes with transforms being in the
///        center of each volume. Also, there's limited capability to transcript
///        boolean volume (Subtraction/Union/Shift). The set of volumes to
///        convert can be constrained using nameList of the Configuration
///        object. Additionally, the factory can also convert the shapes
///        associated with FullPhysicalVolumes convert into bounding volumes.
///        All sufraces constructed from such a volume are then also put into
///        the volume. For compatibility reasons, the conversion constructs
///        Volumes -> processed into TrackingVolume by the Gen1/ Gen 3
///        TrackingGeometry builder and also detector volumes to convert them
///        into Gen2 tracking geometry volumes
class GeoModelDetectorObjectFactory {
 public:
  /// @brief abrivation of the smart pointer to a full physical volume
  using FpvConstLink = GeoModelTree::FpvConstLink;
  ///  @brief Tuple describing the shared ptr to a Volume which will be turned into a TrackingVolume,
  ///          a Gen-2 volume and the pointer to the full physical volume
  using GeoModelVolumeFPVTuple =
      std::tuple<std::shared_ptr<Volume>,
                 std::shared_ptr<Experimental::DetectorVolume>, FpvConstLink>;

  struct Options {
    std::vector<std::string> queries = {};
  };

  // substring matching for additional parameters
  // empty vectors will be matched to everything
  struct Config {
    // List for names to match
    std::vector<std::string> nameList{};
    // List for materials to match
    std::vector<std::string> materialList{};

    // boolean flag to build subvolumes
    bool convertSubVolumes = false;

    // flag to build the desired bounding Boxes
    std::vector<std::string> convertBox;
  };

  struct Cache {
    // The created detector elements and their surfaces
    std::vector<GeoModelSensitiveSurface> sensitiveSurfaces;
    /// @brief Pointer to the surface bound factory
    std::shared_ptr<SurfaceBoundFactory> surfBoundFactory =
        std::make_shared<SurfaceBoundFactory>();
    /// @brief Pointer to the volume bound factory */
    std::shared_ptr<VolumeBoundFactory> volumeBoundFactory =
        std::make_shared<VolumeBoundFactory>();

    // The created representation of bounding boxes  and the corresponding Full
    // Physical Volumes
    std::vector<GeoModelVolumeFPVTuple> volumeBoxFPVs{};
  };

  explicit GeoModelDetectorObjectFactory(
      const Config& cfg,
      std::unique_ptr<const Logger> mlogger = getDefaultLogger(
          "GeoModelDetectorObjectFactory", Acts::Logging::WARNING));

  /// @brief Run the translation from the GeoModelTree to the (sensitive) surfaces.
  /// @param cache: Cache object which will contain the surfaces & box volume bounds
  /// @param gctx: Instance to an geometry context in order to build the envelope volumes
  /// @param geoModelTree: Configured instance of the GeoModelTree to run the construction on
  /// @param options: Options configuring which volumes / materials shall be converted to surfaces
  void construct(Cache& cache, const GeometryContext& gctx,
                 const GeoModelTree& geoModelTree, const Options& options);

  /// @brief Convert a full physical volume (and the appropriate children) into sensitive surfaces
  /// @param name: Published name of the full physical volume in the GeoModelTree
  /// @param fpv: Pointer to the full physical volume to convert
  /// @param cache: Output cache object in which the constructed surfaces are saved
  /// @param gctx: Instance to an geometry context in order to build the envelope volumes
  void convertFpv(const std::string& name, const FpvConstLink& fpv,
                  Cache& cache, const GeometryContext& gctx);

 private:
  ///  @brief Convert the GeoPhysVol into a sensitive Acts::Surface.
  ///  @param geoPV: Pointer to the GeoPhysVol to convert
  ///  @param transform: Placement of the resulting surface in the world
  ///  @param boundFactory: Reference to the BoundFactory to avoid duplicated bounds
  ///                       across similar surfaces
  ///  @param sensitives: Output vector into which the new converted surface is pushed
  void convertSensitive(const PVConstLink& geoPV,
                        const Acts::Transform3& transform,
                        SurfaceBoundFactory& boundFactory,
                        std::vector<GeoModelSensitiveSurface>& sensitives);
  /// @brief Find all sub volumes of a passed volume that are
  ///        good for sensitive detector conversion
  /// @param vol: Pointer to the GeoPhysVol to search through
  /// @return A vector of GeoChildNodeWithTrf containing the information about the
  ///         volumes to convert and their placement w.r.t. the passed volume
  std::vector<GeoChildNodeWithTrf> findAllSubVolumes(
      const PVConstLink& vol) const;
  ///  @brief Checks whether the volume name satisfies the user-defined tokens and/or
  ///         the material of physical volume does it.
  ///  @param name: Name of the physical volume to test. Usually, it's the GeoNameTag or
  ///               the published GeoFullPhysVol entry
  ///  @param physvol: Reference to the PhysicalVolume to additionally check material compatibility
  bool matches(const std::string& name, const PVConstLink& physvol) const;
  ///  @brief Returns whether the name of the published full physical volume is on the list
  ///         to also convert the volume to an envelope volume
  ///  @param name: Name of the published full physical volume
  bool convertBox(const std::string& name) const;

  std::unique_ptr<const Logger> m_logger;
  Config m_cfg;

  const Logger& logger() const { return *m_logger; }
};
}  // namespace Acts

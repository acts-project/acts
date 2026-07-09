// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Utilities/BoundFactory.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/GeoModel/GeoModelDetectorElement.hpp"
#include "ActsPlugins/GeoModel/GeoModelToDetectorVolume.hpp"
#include "ActsPlugins/GeoModel/GeoModelTree.hpp"

#include <GeoModelHelpers/getChildNodesWithTrf.h>

#include "GeoModelKernel/GeoDefinitions.h"

class GeoShape;

namespace Acts {
class Surface;
}  // namespace Acts

namespace ActsPlugins {

/// @addtogroup geomodel_plugin
/// @{

/// Factory class to convert GeoModel objects into Acts volumes and surfaces.
/// The surface conversion
///        process constructs surfaces from GeoTrd, GeoBox, GeoTube,
///        GeoSimplePolygonBrep volume shapes with transforms being in the
///        center of each volume. Also, there's limited capability to transcript
///        boolean volume (Subtraction/Union/Shift). The set of volumes to
///        convert can be constrained using nameList of the Configuration
///        object. Additionally, the factory can also convert the shapes
///        associated with FullPhysicalVolumes into bounding volumes.
///        All sufraces constructed from such a volume are then also put into
///        the volume.
class GeoModelDetectorObjectFactory {
 public:
  /// abrivation of the smart pointer to a full physical volume
  using FpvConstLink = GeoModelTree::FpvConstLink;
  /// Struct containing the link to the full physical volume
  /// that is converted into a tracking volume, the tracking
  /// volume and the associated surfaces.
  struct ConvertedGeoVol {
    /// Pointer to the tracking volume built
    std::unique_ptr<Acts::TrackingVolume> volume{};
    /// Pointer to the full physical volume from which the
    ///        Volumes have been built
    FpvConstLink fullPhysVol{};
    /// List of surfaces belonging to this volume
    std::vector<std::shared_ptr<Acts::Surface>> surfaces{};
  };
  /// List of converted GeoModel volumes
  using ConvertedVolList_t = std::vector<ConvertedGeoVol>;

  /// Options for the GeoModel detector object factory
  struct Options {
    /// Query strings for volume selection
    std::vector<std::string> queries{};
  };

  // substring matching for additional parameters
  // empty vectors will be matched to everything
  /// Configuration for GeoModel detector object factory
  struct Config {
    /// List of name patterns to match
    std::vector<std::string> nameList{};
    /// List of material patterns to match
    std::vector<std::string> materialList{};
    /// Boolean flag to build subvolumes
    bool convertSubVolumes = false;
    /// Flag to build the desired bounding boxes
    std::vector<std::string> convertBox{};
  };

  /// Cache for storing created detector elements and factory objects
  struct Cache {
    /// Pointer to the surface bound factory
    std::shared_ptr<Acts::SurfaceBoundFactory> surfBoundFactory =
        std::make_shared<Acts::SurfaceBoundFactory>();
    /// Pointer to the volume bound factory */
    std::shared_ptr<Acts::VolumeBoundFactory> volumeBoundFactory =
        std::make_shared<Acts::VolumeBoundFactory>();
    /// The created representation of bounding boxes and the corresponding Full
    /// Physical Volumes
    ConvertedVolList_t volumeBoxFPVs{};
  };

  /// Constructor
  /// @param cfg Configuration
  /// @param mlogger Logger
  explicit GeoModelDetectorObjectFactory(
      const Config& cfg,
      std::unique_ptr<const Acts::Logger> mlogger = Acts::getDefaultLogger(
          "GeoModelDetectorObjectFactory", Acts::Logging::WARNING));

  /// Run the translation from the GeoModelTree to the (sensitive) surfaces.
  /// @param cache: Cache object which will contain the surfaces & box volume bounds
  /// @param gctx: Instance to an geometry context in order to build the envelope volumes
  /// @param geoModelTree: Configured instance of the GeoModelTree to run the construction on
  /// @param options: Options configuring which volumes / materials shall be converted to surfaces
  void construct(Cache& cache, const Acts::GeometryContext& gctx,
                 const GeoModelTree& geoModelTree, const Options& options);

  /// Convert a full physical volume (and the appropriate children) into
  /// sensitive surfaces
  /// @param name: Published name of the full physical volume in the GeoModelTree
  /// @param fpv: Pointer to the full physical volume to convert
  /// @param cache: Output cache object in which the constructed surfaces are saved
  /// @param gctx: Instance to an geometry context in order to build the envelope volumes
  void convertFpv(const std::string& name, const FpvConstLink& fpv,
                  Cache& cache, const Acts::GeometryContext& gctx);

 private:
  /// Convert the GeoPhysVol into a sensitive Acts::Surface.
  ///  @param geoPV: Pointer to the GeoPhysVol to convert
  ///  @param transform: Placement of the resulting surface in the world
  ///  @param boundFactory: Reference to the BoundFactory to avoid duplicated bounds
  ///                       across similar surfaces
  ///  @param sensitives: Output vector into which the new converted surface is pushed
  void convertSensitive(const PVConstLink& geoPV,
                        const Acts::Transform3& transform,
                        Acts::SurfaceBoundFactory& boundFactory,
                        std::vector<GeoModelSensitiveSurface>& sensitives);
  /// Find all sub volumes of a passed volume that are
  ///        good for sensitive detector conversion
  /// @param vol: Pointer to the GeoPhysVol to search through
  /// @return A vector of GeoChildNodeWithTrf containing the information about the
  ///         volumes to convert and their placement w.r.t. the passed volume
  std::vector<GeoChildNodeWithTrf> findAllSubVolumes(
      const PVConstLink& vol) const;
  /// Checks whether the volume name satisfies the user-defined tokens and/or
  ///         the material of physical volume does it.
  ///  @param name: Name of the physical volume to test. Usually, it's the GeoNameTag or
  ///               the published GeoFullPhysVol entry
  ///  @param physvol: Reference to the PhysicalVolume to additionally check material compatibility
  bool matches(const std::string& name, const PVConstLink& physvol) const;
  /// Returns whether the name of the published full physical volume is on the
  /// list
  ///         to also convert the volume to an envelope volume
  ///  @param name: Name of the published full physical volume
  bool convertBox(const std::string& name) const;

  std::unique_ptr<const Acts::Logger> m_logger;
  Config m_cfg;

  const Acts::Logger& logger() const { return *m_logger; }
};

/// @}

}  // namespace ActsPlugins

// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Frustum.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Ray.hpp"

#include <functional>
#include <string>
#include <unordered_map>

#include <boost/container/small_vector.hpp>

namespace Acts {

class GlueVolumesDescriptor;
class VolumeBounds;

template <typename object_t>
struct NavigationOptions;

// master typedefs
using TrackingVolumePtr = std::shared_ptr<const TrackingVolume>;
using MutableTrackingVolumePtr = std::shared_ptr<TrackingVolume>;

using TrackingVolumeBoundaryPtr =
    std::shared_ptr<const BoundarySurfaceT<TrackingVolume>>;
using TrackingVolumeBoundaries = std::vector<TrackingVolumeBoundaryPtr>;

// possible contained
using TrackingVolumeArray = BinnedArray<TrackingVolumePtr>;
using TrackingVolumeVector = std::vector<TrackingVolumePtr>;
using MutableTrackingVolumeVector = std::vector<MutableTrackingVolumePtr>;
using LayerArray = BinnedArray<LayerPtr>;
using LayerVector = std::vector<LayerPtr>;

// Intersection with Layer
using LayerIntersection = ObjectIntersection<Layer, Surface>;

/// BoundarySurface of a volume
using BoundarySurface = BoundarySurfaceT<TrackingVolume>;
/// Intersection with a @c BoundarySurface
using BoundaryIntersection = ObjectIntersection<BoundarySurface, Surface>;

/// @class TrackingVolume
///
/// Full Volume description used in Tracking,
/// it inherits from Volume to get the geometrical structure.
///
///     A TrackingVolume at navigation level can provide the (layer) material
/// information / internal navigation with in
///     5 different ways:
///
///         --- a) Static confinement of Layers
///         --- b) detached sub volumes
///         --- b) unordered (arbitrarily oriented) layers
///         --- d) unordered sub volumes
///         --- e) unordered layers AND unordered subvolumes
///
///    The TrackingVolume can also be a simple container of other
///    TrackingVolumes
///
/// In addition it is capable of holding a subarray of Layers and
/// TrackingVolumes.
///
class TrackingVolume : public Volume {
  friend class TrackingGeometry;

 public:
  TrackingVolume() = delete;
  ~TrackingVolume() override;
  TrackingVolume(const TrackingVolume&) = delete;
  TrackingVolume& operator=(const TrackingVolume&) = delete;

  /// Factory constructor for a container TrackingVolume
  /// - by definition a Vacuum volume
  ///
  /// @param transform is the global 3D transform to position the volume in
  /// space
  /// @param volumeBounds is the description of the volume boundaries
  /// @param containedVolumes are the static volumes that fill this volume
  /// @param volumeName is a string identifier
  ///
  /// @return shared pointer to a new TrackingVolume
  static MutableTrackingVolumePtr create(
      const Transform3& transform, VolumeBoundsPtr volumeBounds,
      const std::shared_ptr<const TrackingVolumeArray>& containedVolumes =
          nullptr,
      const std::string& volumeName = "undefined") {
    return MutableTrackingVolumePtr(new TrackingVolume(
        transform, std::move(volumeBounds), containedVolumes, volumeName));
  }

  /// Factory constructor for Tracking Volume with a bounding volume hierarchy
  ///
  /// @param transform is the global 3D transform to position the volume in
  /// space
  /// @param volbounds is the description of the volume boundaries
  /// @param boxStore Vector owning the contained bounding boxes
  /// @param descendants Vector owning the child volumes
  /// @param top The top of the hierarchy (top node)
  /// @param volumeMaterial is the materials of the tracking volume
  /// @param volumeName is a string identifier
  ///
  /// @return shared pointer to a new TrackingVolume
  static MutableTrackingVolumePtr create(
      const Transform3& transform, VolumeBoundsPtr volbounds,
      std::vector<std::unique_ptr<Volume::BoundingBox>> boxStore,
      std::vector<std::unique_ptr<const Volume>> descendants,
      const Volume::BoundingBox* top,
      std::shared_ptr<const IVolumeMaterial> volumeMaterial,
      const std::string& volumeName = "undefined") {
    return MutableTrackingVolumePtr(new TrackingVolume(
        transform, std::move(volbounds), std::move(boxStore),
        std::move(descendants), top, std::move(volumeMaterial), volumeName));
  }

  /// Factory constructor for Tracking Volumes with content
  /// - can not be a container volume
  ///
  /// @param transform is the global 3D transform to position the volume in
  /// space
  /// @param volumeBounds is the description of the volume boundaries
  /// @param volumeMaterial is are materials of the tracking volume
  /// @param containedLayers is the confined layer array (optional)
  /// @param containedVolumes is the confined volume array (optional)
  /// @param denseVolumes is the array of dense volulmes (optional)
  /// @param volumeName is a string identifier
  ///
  /// @return shared pointer to a new TrackingVolume
  static MutableTrackingVolumePtr create(
      const Transform3& transform, VolumeBoundsPtr volumeBounds,
      std::shared_ptr<const IVolumeMaterial> volumeMaterial,
      std::unique_ptr<const LayerArray> containedLayers = nullptr,
      std::shared_ptr<const TrackingVolumeArray> containedVolumes = nullptr,
      MutableTrackingVolumeVector denseVolumes = {},
      const std::string& volumeName = "undefined") {
    return MutableTrackingVolumePtr(new TrackingVolume(
        transform, std::move(volumeBounds), std::move(volumeMaterial),
        std::move(containedLayers), std::move(containedVolumes),
        std::move(denseVolumes), volumeName));
  }

  /// Return the associated Layer to the global position
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the associated global position
  ///
  /// @return plain pointer to layer object
  const Layer* associatedLayer(const GeometryContext& gctx,
                               const Vector3& position) const;

  /// @brief Resolves the volume into (compatible) Layers
  ///
  /// This is the method for the propagator/extrapolator
  /// @tparam options_t Type of navigation options object for decomposition
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position Position for the search
  /// @param direction Direction for the search
  /// @param options The templated navigation options
  ///
  /// @return vector of compatible intersections with layers
  boost::container::small_vector<LayerIntersection, 10> compatibleLayers(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction, const NavigationOptions<Layer>& options) const;

  /// @brief Returns all boundary surfaces sorted by the user.
  ///
  /// @tparam options_t Type of navigation options object for decomposition
  /// @tparam sorter_t Type of the boundary surface sorter
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position for searching
  /// @param direction The direction for searching
  /// @param options The templated navigation options
  /// @param logger A @c LoggerWrapper instance
  ///
  /// @return is the templated boundary intersection
  boost::container::small_vector<BoundaryIntersection, 4> compatibleBoundaries(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction, const NavigationOptions<Surface>& options,
      LoggerWrapper logger = getDummyLogger()) const;

  /// @brief Return surfaces in given direction from bounding volume hierarchy
  /// @tparam options_t Type of navigation options object for decomposition
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position to start from
  /// @param direction The direction towards which to test
  /// @param angle The opening angle
  /// @param options The templated navigation options
  ///
  /// @return Vector of surface candidates
  std::vector<SurfaceIntersection> compatibleSurfacesFromHierarchy(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction, double angle,
      const NavigationOptions<Surface>& options) const;

  /// Return the associated sub Volume, returns THIS if no subVolume exists
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global position associated with that search
  /// @param tol Search position tolerance for dense volumes
  ///
  /// @return plain pointer to associated with the position
  const TrackingVolume* lowestTrackingVolume(const GeometryContext& gctx,
                                             const Vector3& position,
                                             const double tol = 0.) const;

  /// Return the confined static layer array - if it exists
  /// @return the BinnedArray of static layers if exists
  const LayerArray* confinedLayers() const;

  /// Return the confined volumes of this container array - if it exists
  std::shared_ptr<const TrackingVolumeArray> confinedVolumes() const;

  /// Return the confined dense volumes
  const MutableTrackingVolumeVector denseVolumes() const;

  /// @brief Visit all sensitive surfaces
  ///
  /// @tparam visitor_t Type of the callable visitor
  ///
  /// @param visitor The callable. Will be called for each sensitive surface
  /// that is found
  ///
  /// If a context is needed for the vist, the vistitor has to provide this
  /// e.g. as a private member
  template <typename visitor_t>
  void visitSurfaces(visitor_t&& visitor) const {
    if (!m_confinedVolumes) {
      // no sub volumes => loop over the confined layers
      if (m_confinedLayers) {
        for (const auto& layer : m_confinedLayers->arrayObjects()) {
          if (layer->surfaceArray() == nullptr) {
            // no surface array (?)
            continue;
          }
          for (const auto& srf : layer->surfaceArray()->surfaces()) {
            visitor(srf);
          }
        }
      }
    } else {
      // contains sub volumes
      for (const auto& volume : m_confinedVolumes->arrayObjects()) {
        volume->visitSurfaces(visitor);
      }
    }
  }

  /// Returns the VolumeName - for debug reason, might be depreciated later
  const std::string& volumeName() const;

  /// Method to return the BoundarySurfaces
  const TrackingVolumeBoundaries& boundarySurfaces() const;

  /// Return the material of the volume
  const IVolumeMaterial* volumeMaterial() const;

  /// Return the material of the volume as shared pointer
  const std::shared_ptr<const IVolumeMaterial>& volumeMaterialSharedPtr() const;

  /// Set the boundary surface material description
  ///
  /// The material is usually derived in a complicated way and loaded from
  /// a framework given source. As various volumes could potentially share the
  /// the same material description, it is provided as a shared object
  ///
  /// @param surfaceMaterial Material description of this volume
  /// @param bsFace Specifies which boundary surface to assign the material to
  void assignBoundaryMaterial(
      std::shared_ptr<const ISurfaceMaterial> surfaceMaterial,
      BoundarySurfaceFace bsFace);

  /// Set the volume material description
  ///
  /// The material is usually derived in a complicated way and loaded from
  /// a framework given source. As various volumes could potentially share the
  /// the same material description, it is provided as a shared object
  ///
  /// @param material Material description of this volume
  void assignVolumeMaterial(std::shared_ptr<const IVolumeMaterial> material);

  /// Glue another tracking volume to this one
  ///  - if common face is set the glued volumes are sharing the boundary, down
  /// to the last navigation volume
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bsfMine is the boundary face indicater where to glue
  /// @param neighbor is the TrackingVolume to be glued
  /// @param bsfNeighbor is the boudnary surface of the neighbor
  void glueTrackingVolume(const GeometryContext& gctx,
                          BoundarySurfaceFace bsfMine, TrackingVolume* neighbor,
                          BoundarySurfaceFace bsfNeighbor);

  /// Glue another tracking volume to this one
  ///  - if common face is set the glued volumes are sharing the boundary, down
  /// to the last navigation volume
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bsfMine is the boundary face indicater where to glue
  /// @param neighbors are the TrackingVolumes to be glued
  /// @param bsfNeighbor are the boudnary surface of the neighbors
  void glueTrackingVolumes(
      const GeometryContext& gctx, BoundarySurfaceFace bsfMine,
      const std::shared_ptr<TrackingVolumeArray>& neighbors,
      BoundarySurfaceFace bsfNeighbor);

  /// Provide a new BoundarySurface from the glueing
  ///
  /// @param bsf is the boundary face indicater where to glue
  /// @param bs is the new boudnary surface
  /// @param checkmaterial is a flag how to deal with material, if true:
  /// - if the old boundary surface had a material description
  ///   but the new one has not, keep the current one
  /// - in all other cases just assign the new boundary surface
  void updateBoundarySurface(
      BoundarySurfaceFace bsf,
      std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bs,
      bool checkmaterial = true);

  /// Register the outside glue volumes -
  /// ordering is in the TrackingVolume Frame:
  ///  - negativeFaceXY
  ///  - (faces YZ, ZY, radial faces)
  ///  - positiveFaceXY
  ///
  /// @param gvd register a new GlueVolumeDescriptor
  /// @todo update to shared/unique ptr
  void registerGlueVolumeDescriptor(GlueVolumesDescriptor* gvd);

  /// Register the outside glue volumes -
  /// ordering is in the TrackingVolume Frame:
  ///  - negativeFaceXY
  ///  - (faces YZ, ZY, radial faces)
  ///  - positiveFaceXY
  GlueVolumesDescriptor& glueVolumesDescriptor();

  /// Return whether this TrackingVolume has a BoundingVolumeHierarchy
  /// associated
  /// @return If it has a BVH or not.
  bool hasBoundingVolumeHierarchy() const;

  /// Register the color code
  ///
  /// @param icolor is a color number
  void registerColorCode(unsigned int icolor);

  /// Get the color code
  unsigned int colorCode() const;

  /// Return the MotherVolume - if it exists
  const TrackingVolume* motherVolume() const;

  /// Set the MotherVolume
  ///
  /// @param mvol is the mother volume
  void setMotherVolume(const TrackingVolume* mvol);

 protected:
  /// Constructor for a container Volume
  /// - vacuum filled volume either as a for other tracking volumes
  ///
  /// @param transform is the global 3D transform to position the volume in
  /// space
  /// @param volbounds is the description of the volume boundaries
  /// @param containedVolumeArray are the static volumes that fill this volume
  /// @param volumeName is a string identifier
  TrackingVolume(const Transform3& transform, VolumeBoundsPtr volbounds,
                 const std::shared_ptr<const TrackingVolumeArray>&
                     containedVolumeArray = nullptr,
                 const std::string& volumeName = "undefined");

  TrackingVolume(const Transform3& transform, VolumeBoundsPtr volbounds,
                 std::vector<std::unique_ptr<Volume::BoundingBox>> boxStore,
                 std::vector<std::unique_ptr<const Volume>> descendants,
                 const Volume::BoundingBox* top,
                 std::shared_ptr<const IVolumeMaterial> volumeMaterial,
                 const std::string& volumeName = "undefined");

  /// Constructor for a full equipped Tracking Volume
  ///
  /// @param transform is the global 3D transform to position the volume in
  /// space
  /// @param volumeBounds is the description of the volume boundaries
  /// @param volumeMaterial is are materials of the tracking volume
  /// @param staticLayerArray is the confined layer array (optional)
  /// @param containedVolumeArray are the sub volumes if the volume is a
  /// container
  /// @param denseVolumeVector  The contained dense volumes
  /// @param volumeName is a string identifier
  TrackingVolume(
      const Transform3& transform, VolumeBoundsPtr volumeBounds,
      std::shared_ptr<const IVolumeMaterial> volumeMaterial,
      std::unique_ptr<const LayerArray> staticLayerArray = nullptr,
      std::shared_ptr<const TrackingVolumeArray> containedVolumeArray = nullptr,
      MutableTrackingVolumeVector denseVolumeVector = {},
      const std::string& volumeName = "undefined");

 private:
  void connectDenseBoundarySurfaces(
      MutableTrackingVolumeVector& confinedDenseVolumes);

  /// Create Boundary Surface
  void createBoundarySurfaces();

  /// method to synchronize the layers with potentially updated volume bounds:
  /// - adapts the layer dimensions to the new volumebounds + envelope
  ///
  /// @param envelope is the clearance between volume boundary and layer
  void synchronizeLayers(double envelope = 1.) const;

  /// close the Geometry, i.e. set the GeometryIdentifier and assign material
  ///
  /// @param materialDecorator is a dedicated decorator for the
  ///        material to be assigned (surface, volume based)
  /// @param volumeMap is a map to find the a volume by identifier
  /// @param vol is the geometry id of the volume
  ///        as calculated by the TrackingGeometry
  /// @param hook Identifier hook to be applied to surfaces
  ///
  void closeGeometry(
      const IMaterialDecorator* materialDecorator,
      std::unordered_map<GeometryIdentifier, const TrackingVolume*>& volumeMap,
      size_t& vol, const GeometryIdentifierHook& hook);

  /// interlink the layers in this TrackingVolume
  void interlinkLayers();

  /// The volume based material the TrackingVolume consists of
  std::shared_ptr<const IVolumeMaterial> m_volumeMaterial{nullptr};

  /// Remember the mother volume
  const TrackingVolume* m_motherVolume{nullptr};

  // the boundary surfaces
  std::vector<TrackingVolumeBoundaryPtr> m_boundarySurfaces;

  ///(a) static configuration ordered by Binned arrays
  /// static layers
  std::unique_ptr<const LayerArray> m_confinedLayers = nullptr;

  /// Array of Volumes inside the Volume when actin as container
  std::shared_ptr<const TrackingVolumeArray> m_confinedVolumes = nullptr;

  /// confined dense
  MutableTrackingVolumeVector m_confinedDenseVolumes;

  /// Volumes to glue Volumes from the outside
  GlueVolumesDescriptor* m_glueVolumeDescriptor{nullptr};

  /// Volume name for debug reasons & screen output
  std::string m_name;

  /// color code for displaying
  unsigned int m_colorCode{20};

  /// Bounding Volume Hierarchy (BVH)
  std::vector<std::unique_ptr<const Volume::BoundingBox>> m_boundingBoxes;
  std::vector<std::unique_ptr<const Volume>> m_descendantVolumes;
  const Volume::BoundingBox* m_bvhTop{nullptr};
};

inline const std::string& TrackingVolume::volumeName() const {
  return m_name;
}

inline const IVolumeMaterial* TrackingVolume::volumeMaterial() const {
  return m_volumeMaterial.get();
}

inline const std::shared_ptr<const IVolumeMaterial>&
TrackingVolume::volumeMaterialSharedPtr() const {
  return m_volumeMaterial;
}

inline void TrackingVolume::assignVolumeMaterial(
    std::shared_ptr<const IVolumeMaterial> material) {
  m_volumeMaterial = std::move(material);
}

inline const LayerArray* TrackingVolume::confinedLayers() const {
  return m_confinedLayers.get();
}

inline const MutableTrackingVolumeVector TrackingVolume::denseVolumes() const {
  return m_confinedDenseVolumes;
}

inline std::shared_ptr<const TrackingVolumeArray>
TrackingVolume::confinedVolumes() const {
  return m_confinedVolumes;
}

inline void TrackingVolume::registerColorCode(unsigned int icolor) {
  m_colorCode = icolor;
}

inline unsigned int TrackingVolume::colorCode() const {
  return m_colorCode;
}

inline const TrackingVolume* TrackingVolume::motherVolume() const {
  return m_motherVolume;
}

inline void TrackingVolume::setMotherVolume(const TrackingVolume* mvol) {
  m_motherVolume = mvol;
}

inline bool TrackingVolume::hasBoundingVolumeHierarchy() const {
  return m_bvhTop != nullptr;
}

#ifndef DOXYGEN
#include "Acts/Geometry/detail/TrackingVolume.ipp"
#endif

}  // namespace Acts

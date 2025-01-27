// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/BoundarySurfaceFace.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingVolumeVisitorConcept.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Surfaces/SurfaceVisitorConcept.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Concepts.hpp"
#include "Acts/Utilities/Frustum.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Ray.hpp"

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts {

class GlueVolumesDescriptor;
class VolumeBounds;
template <typename object_t>
struct NavigationOptions;
class GeometryIdentifier;
class IMaterialDecorator;
class ISurfaceMaterial;
class IVolumeMaterial;
class Surface;
class TrackingVolume;
struct GeometryIdentifierHook;

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

/// Intersection with @c Layer
using LayerIntersection = std::pair<SurfaceIntersection, const Layer*>;
/// Multi-intersection with @c Layer
using LayerMultiIntersection =
    std::pair<SurfaceMultiIntersection, const Layer*>;

/// BoundarySurface of a volume
using BoundarySurface = BoundarySurfaceT<TrackingVolume>;
/// Intersection with a @c BoundarySurface
using BoundaryIntersection =
    std::pair<SurfaceIntersection, const BoundarySurface*>;
/// Multi-intersection with a @c BoundarySurface
using BoundaryMultiIntersection =
    std::pair<SurfaceMultiIntersection, const BoundarySurface*>;

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

  /// Constructor for a container Volume
  /// - vacuum filled volume either as a for other tracking volumes
  ///
  /// @param transform is the global 3D transform to position the volume in
  /// space
  /// @param volbounds is the description of the volume boundaries
  /// @param containedVolumeArray are the static volumes that fill this volume
  /// @param volumeName is a string identifier
  TrackingVolume(const Transform3& transform,
                 std::shared_ptr<VolumeBounds> volbounds,
                 const std::shared_ptr<const TrackingVolumeArray>&
                     containedVolumeArray = nullptr,
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
      const Transform3& transform, std::shared_ptr<VolumeBounds> volumeBounds,
      std::shared_ptr<const IVolumeMaterial> volumeMaterial,
      std::unique_ptr<const LayerArray> staticLayerArray = nullptr,
      std::shared_ptr<const TrackingVolumeArray> containedVolumeArray = nullptr,
      MutableTrackingVolumeVector denseVolumeVector = {},
      const std::string& volumeName = "undefined");

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
  /// @param logger A @c Logger instance
  ///
  /// @return is the templated boundary intersection
  boost::container::small_vector<BoundaryIntersection, 4> compatibleBoundaries(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction, const NavigationOptions<Surface>& options,
      const Logger& logger = getDummyLogger()) const;

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

  /// @brief Visit all reachable surfaces
  ///
  /// @tparam visitor_t Type of the callable visitor
  ///
  /// @param visitor The callable. Will be called for each reachable surface
  /// that is found, a selection of the surfaces can be done in the visitor
  /// @param restrictToSensitives If true, only sensitive surfaces are visited
  ///
  /// @note If a context is needed for the visit, the vistitor has to provide
  /// this, e.g. as a private member
  template <ACTS_CONCEPT(SurfaceVisitor) visitor_t>
  void visitSurfaces(visitor_t&& visitor, bool restrictToSensitives) const {
    if (!restrictToSensitives) {
      // Visit the boundary surfaces
      for (const auto& bs : m_boundarySurfaces) {
        visitor(&(bs->surfaceRepresentation()));
      }
    }

    // Internal structure
    if (m_confinedVolumes == nullptr) {
      // no sub volumes => loop over the confined layers
      if (m_confinedLayers != nullptr) {
        for (const auto& layer : m_confinedLayers->arrayObjects()) {
          // Surfaces contained in the surface array
          if (layer->surfaceArray() != nullptr) {
            for (const auto& srf : layer->surfaceArray()->surfaces()) {
              visitor(srf);
              continue;
            }
          }
          if (!restrictToSensitives) {
            // Surfaces of the layer
            visitor(&layer->surfaceRepresentation());
            // Approach surfaces of the layer
            if (layer->approachDescriptor() != nullptr) {
              for (const auto& srf :
                   layer->approachDescriptor()->containedSurfaces()) {
                visitor(srf);
              }
            }
          }
        }
      }
    } else {
      // contains sub volumes
      for (const auto& volume : m_confinedVolumes->arrayObjects()) {
        volume->visitSurfaces(visitor, restrictToSensitives);
      }
    }
  }

  /// @brief Visit all sensitive surfaces
  ///
  /// @tparam visitor_t Type of the callable visitor
  ///
  /// @param visitor The callable. Will be called for each sensitive surface
  /// that is found, a selection of the surfaces can be done in the visitor
  ///
  /// @note If a context is needed for the visit, the vistitor has to provide
  /// this, e.g. as a private member
  template <ACTS_CONCEPT(SurfaceVisitor) visitor_t>
  void visitSurfaces(visitor_t&& visitor) const {
    visitSurfaces(std::forward<visitor_t>(visitor), true);
  }

  /// @brief Visit all reachable tracking volumes
  ///
  /// @tparam visitor_t Type of the callable visitor
  ///
  /// @param visitor The callable. Will be called for each reachable volume
  /// that is found, a selection of the volumes can be done in the visitor
  ///
  /// @note If a context is needed for the visit, the vistitor has to provide
  /// this, e.g. as a private member
  template <ACTS_CONCEPT(TrackingVolumeVisitor) visitor_t>
  void visitVolumes(visitor_t&& visitor) const {
    visitor(this);
    if (m_confinedVolumes != nullptr) {
      // contains sub volumes
      for (const auto& volume : m_confinedVolumes->arrayObjects()) {
        volume->visitVolumes(visitor);
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
  /// @param bsfNeighbor is the boundary surface of the neighbor
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
  /// @param bsfNeighbor are the boundary surface of the neighbors
  void glueTrackingVolumes(
      const GeometryContext& gctx, BoundarySurfaceFace bsfMine,
      const std::shared_ptr<TrackingVolumeArray>& neighbors,
      BoundarySurfaceFace bsfNeighbor);

  /// Provide a new BoundarySurface from the glueing
  ///
  /// @param bsf is the boundary face indicater where to glue
  /// @param bs is the new boundary surface
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

  /// Return the MotherVolume - if it exists
  const TrackingVolume* motherVolume() const;

  /// Return the MotherVolume - if it exists
  TrackingVolume* motherVolume();

  /// Set the MotherVolume
  ///
  /// @param mvol is the mother volume
  void setMotherVolume(TrackingVolume* mvol);

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
  /// @param logger A @c Logger instance
  ///
  void closeGeometry(
      const IMaterialDecorator* materialDecorator,
      std::unordered_map<GeometryIdentifier, const TrackingVolume*>& volumeMap,
      std::size_t& vol, const GeometryIdentifierHook& hook,
      const Logger& logger = getDummyLogger());

  /// interlink the layers in this TrackingVolume
  void interlinkLayers();

  /// The volume based material the TrackingVolume consists of
  std::shared_ptr<const IVolumeMaterial> m_volumeMaterial{nullptr};

  /// Remember the mother volume
  TrackingVolume* m_motherVolume{nullptr};

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
};

}  // namespace Acts

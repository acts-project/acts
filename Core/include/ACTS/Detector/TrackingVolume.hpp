// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingVolume.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_TRACKINGVOLUME_H
#define ACTS_DETECTOR_TRACKINGVOLUME_H 1

#include <map>
#include <string>
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Surfaces/BoundaryCheck.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/GeometrySignature.hpp"
#include "ACTS/Volumes/BoundarySurfaceT.hpp"
#include "ACTS/Volumes/Volume.hpp"

namespace Acts {

// classes
class TrackingVolume;
class DetachedTrackingVolume;
class GlueVolumesDescriptor;
class VolumeBounds;
class Material;

typedef std::shared_ptr<const BoundarySurfaceT<TrackingVolume>>
    TrackingVolumeBoundaryPtr;

// master typedefs
typedef std::shared_ptr<const TrackingVolume>         TrackingVolumePtr;
typedef std::shared_ptr<const DetachedTrackingVolume> DetachedTrackingVolumePtr;

// possible contained
typedef BinnedArray<TrackingVolumePtr>         TrackingVolumeArray;
typedef std::vector<TrackingVolumePtr>         TrackingVolumeVector;
typedef BinnedArray<LayerPtr>                  LayerArray;
typedef std::vector<LayerPtr>                  LayerVector;
typedef std::vector<DetachedTrackingVolumePtr> DetachedVolumeVector;

template <class T>
using LayerIntersection = FullIntersection<Layer, Surface, T>;
template <class T>
using BoundaryIntersection
    = FullIntersection<BoundarySurfaceT<TrackingVolume>, Surface, T>;

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
class TrackingVolume : public Volume
{
  friend class TrackingGeometry;

public:
  /// Destructor
  ~TrackingVolume();

  /// Factory constructor for a conatiner TrackingVolume
  /// - by definition a Vacuum volume
  ///
  /// @param htrans is the global 3D transform to position the volume ins pace
  /// @param volumeBounds is the description of the volume boundaries
  /// @param containedVolumes are the static volumes that fill this volume
  /// @param volumeName is a string identifier
  ///
  /// @return shared pointer to a new TrackingVolume
  static TrackingVolumePtr
  create(std::shared_ptr<Transform3D>               htrans,
         VolumeBoundsPtr                            volumeBounds,
         std::shared_ptr<const TrackingVolumeArray> containedVolumes = nullptr,
         const std::string&                         volumeName = "undefined")
  {
    return TrackingVolumePtr(
        new TrackingVolume(htrans, volumeBounds, containedVolumes, volumeName));
  }

  /// Factory constructor for Tracking Volumes with content
  /// - can not be a container volume
  ///
  /// @param htrans is the global 3D transform to position the volume ins pace
  /// @param volumeBounds is the description of the volume boundaries
  /// @param matprop is are materials of the tracking volume
  /// @param cLayerArray is the confined layer array (optional)
  /// @param cLayerVector is the confined arbitrary layer vector
  /// @param cVolumeVector is the confined arbitrary volume vector
  /// @param dVolumeVector is the confined arbeitrary detached volume vector
  /// @param volumeName is a string identifier
  ///
  /// @return shared pointer to a new TrackingVolume
  static TrackingVolumePtr
  create(std::shared_ptr<Transform3D>      htrans,
         VolumeBoundsPtr                   volumeBounds,
         std::shared_ptr<Material>         matprop,
         std::unique_ptr<const LayerArray> cLayerArray   = nullptr,
         const LayerVector                 cLayerVector  = {},
         const TrackingVolumeVector        cVolumeVector = {},
         const DetachedVolumeVector        dVolumeVector = {},
         const std::string&                volumeName    = "undefined")
  {
    return TrackingVolumePtr(new TrackingVolume(htrans,
                                                volumeBounds,
                                                matprop,
                                                std::move(cLayerArray),
                                                cLayerVector,
                                                nullptr,
                                                cVolumeVector,
                                                dVolumeVector,
                                                volumeName));
  }
  /// Factory constructor with a shift
  ///
  /// @param trVol is the source tracking volume
  /// @param shift is the additional shift applied after copying
  /// @param volumeName is a string identifier
  ///
  /// @return shared pointer to a new TrackingVolume
  static TrackingVolumePtr
  create(const TrackingVolume& trVol,
         const Transform3D&    shift,
         const std::string&    volumeName = "undefined")
  {
    return TrackingVolumePtr(new TrackingVolume(trVol, shift, volumeName));
  }

  /// Return the associated Layer to the global position
  ///
  /// @param gpos is the associated global position
  ///
  /// @return plain pointer to layer object
  const Layer*
  associatedLayer(const Vector3D& gpos) const;

  /// Return the material layers ordered based on straight line intersections:
  ///
  /// - startLayer and endLayer are included in the list
  /// @param sLayer is the start layer for the search (@todo docu: check if
  /// returned)
  /// @param eLayer is the end layer for the search (@todo docu: check if
  /// returned)
  /// @param parameters are the templated parameters for searching
  /// @param pDir is an additional direction prescription
  /// @param bcheck is a boundary check directive
  /// @param resolveMaterial is the prescription how to deal with material
  /// @param resolveSubSurfaces is the prescription on how to deal with
  /// sensitive surfaces
  ///
  /// @return intersection wiht the layer
  template <class T>
  std::vector<LayerIntersection<T>>
  layerCandidatesOrdered(const Layer*         sLayer,
                         const Layer*         eLayer,
                         const T&             parameters,
                         PropDirection        pDir            = alongMomentum,
                         const BoundaryCheck& bcheck            = true,
                         bool                 resolveMaterial = true,
                         bool                 resolveSubSurfaces = false) const;

  /// Return the associated sub Volume, returns THIS if no subVolume exists
  ///
  /// @param gpos is the global position associated with that search
  ///
  /// @return plain pointer to associated with the position
  const TrackingVolume*
  trackingVolume(const Vector3D& gpos) const;

  /// Return the next volume along the navigation stream
  ///
  /// @param pos is the glboal position associated with that search
  /// @param dir is the direction for this search
  /// @param pDir is an additional direction prescription
  ///
  /// @return the next associated tracking volume, nullptr if it does not exist
  const TrackingVolume*
  nextVolume(const Vector3D& pos,
             const Vector3D& dir,
             PropDirection   pDir = alongMomentum) const;

  /// Return the dynamically created vector of detached sub volumes
  ///
  /// @param pos is the glboal position associated with that search
  /// @param tol is the absolute tolerance for the search           
  ///
  /// @return the list of associated detached tracking volumes, nullptr if it
  /// does not exist
  const DetachedVolumeVector*
  detachedTrackingVolumes(const Vector3D& pos, double tol) const;

  /// Return the confined static layer array - if it exists
  /// @return the BinnedArray of static layers if exists
  const LayerArray*
  confinedLayers() const;

  /// Return the arbitrary layer array - if it exists
  /// @return the vector of arbitrary layers
  const LayerVector
  confinedArbitraryLayers() const;

  /// Return the confined volumes of this container array - if it exists
  std::shared_ptr<const TrackingVolumeArray>
  confinedVolumes() const;

  /// Return the confind volume array as a shared pointer - for glueing
  std::shared_ptr<const TrackingVolumeArray>
  confinedVolumesSharedPtr() const;

  /// Return detached subVolumes - not the ownership
  const DetachedVolumeVector
  confinedDetachedVolumes() const;

  /// Return unordered subVolumes - not the ownership
  const TrackingVolumeVector
  confinedDenseVolumes() const;

  /// Returns the VolumeName - for debug reason, might be depreciated later
  const std::string&
  volumeName() const;

  /// Method to return the BoundarySurfaces
  const std::vector<TrackingVolumeBoundaryPtr>&
  boundarySurfaces() const;

  /// Returns the boundary surfaces ordered in probability to hit them based
  /// on straight line intersection
  ///
  /// @param parameters are the templated tracking parameters
  /// @param pDir is the additional direction presciprion
  ///
  /// @return is the templated boundary intersection
  template <class T>
  std::vector<BoundaryIntersection<T>>
  boundarySurfacesOrdered(const T& parameters,
                          PropDirection pDir = alongMomentum) const;

  /// check if you are on a boundary surface
  ///
  /// @param parameters is the type of track parameters
  ///
  /// @return boolean indicator if the parameters are on boundary
  template <class T>
  bool
  onVolumeBoundary(const T& parameters) const;

  /// Glue another tracking volume to this one
  ///  - if common face is set the glued volumes are sharing the boundary, down
  /// to the last navigation volume
  ///
  /// @param bsfMine is the boundary face indicater where to glue
  /// @param neighbor is the TrackingVolume to be glued
  /// @param bsfNeighbor is the boudnary surface of the neighbor
  void
  glueTrackingVolume(BoundarySurfaceFace bsfMine,
                     TrackingVolumePtr   neighbor,
                     BoundarySurfaceFace bsfNeighbor) const;

  /// Glue another tracking volume to this one
  ///  - if common face is set the glued volumes are sharing the boundary, down
  /// to the last navigation volume
  ///
  ///
  /// @param bsfMine is the boundary face indicater where to glue
  /// @param neighbors are the TrackingVolumes to be glued
  /// @param bsfNeighbors are the boudnary surface of the neighbors
  void
  glueTrackingVolumes(BoundarySurfaceFace                        bsfMine,
                      std::shared_ptr<const TrackingVolumeArray> neighbors,
                      BoundarySurfaceFace bsfNeighbors) const;

  /// Provide a new BoundarySurface from the glueing
  ///
  ///
  /// @param bsfMine is the boundary face indicater where to glue
  /// @param bsSurface is the new boudnary surface
  void
  updateBoundarySurface(
      BoundarySurfaceFace                                     bsfMine,
      std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bsSurface) const;

  /// Register the outside glue volumes -
  /// ordering is in the TrackingVolume Frame:
  ///  - negativeFaceXY
  ///  - (faces YZ, ZY, radial faces)
  ///  - positiveFaceXY
  ///
  /// @param gvd register a new GlueVolumeDescriptor
  /// @todo update to shared/unique ptr
  void
  registerGlueVolumeDescriptor(GlueVolumesDescriptor* gvd) const;

  /// Register the outside glue volumes -
  /// ordering is in the TrackingVolume Frame:
  ///  - negativeFaceXY
  ///  - (faces YZ, ZY, radial faces)
  ///  - positiveFaceXY
  const GlueVolumesDescriptor&
  glueVolumesDescriptor() const;

  /// Sign the volume - the geometry builder has to do that
  ///
  /// @param signat is the volume signature
  /// @param geotype is the volume navigation type
  void
  sign(GeometrySignature signat, GeometryType geotype = Static) const;

  /// return the Signature
  GeometrySignature
  geometrySignature() const;

  /// return the Signature
  GeometryType
  geometryType() const;

  /// Register the color code
  ///
  /// @param icolor is a color number
  void
  registerColorCode(unsigned int icolor) const;

  /// Get the color code
  unsigned int
  colorCode() const;

  /// Return the MotherVolume - if it exists
  const TrackingVolume*
  motherVolume() const;

  /// Return the MotherVolume - if it exists
  ///
  /// @param mvol is the mother volume
  void
  setMotherVolume(const TrackingVolume* mvol) const;

  /// Add Material
  ///
  /// @param mat is the materiall to be added
  /// @param fact is the fraction factor
  void
  addMaterial(std::shared_ptr<const Material> mat, float fact = 1.);

protected:
  /// Default constructor
  TrackingVolume();

  /// Constructor for a container Volume
  /// - vacuum filled volume either as a for other tracking volumes
  ///
  /// @param htrans is the global 3D transform to position the volume ins pace
  /// @param volumeBounds is the description of the volume boundaries
  /// @param containedVolumes are the static volumes that fill this volume
  /// @param volumeName is a string identifier
  TrackingVolume(
      std::shared_ptr<Transform3D>                     htrans,
      VolumeBoundsPtr                                  volumeBounds,
      const std::shared_ptr<const TrackingVolumeArray> containedVolumes
      = nullptr,
      const std::string& volumeName = "undefined");

  /// Constructor for a full equipped Tracking Volume
  ///
  /// @param htrans is the global 3D transform to position the volume ins pace
  /// @param volumeBounds is the description of the volume boundaries
  /// @param matprop is are materials of the tracking volume
  /// @param cLayerArray is the confined layer array (optional)
  /// @param cLayerVector is the confined arbitrary layer vector
  /// @param cVolumeArray is the confined volume array
  /// @param cVolumeVector is the confined arbitrary volume vector
  /// @param dVolumeVector is the confined arbeitrary detached volume vector
  /// @param volumeName is a string identifier
  TrackingVolume(std::shared_ptr<Transform3D>      htrans,
                 VolumeBoundsPtr                   volumeBounds,
                 std::shared_ptr<Material>         matprop,
                 std::unique_ptr<const LayerArray> cLayerArray  = nullptr,
                 const LayerVector                 cLayerVector = {},
                 std::shared_ptr<const TrackingVolumeArray> cVolumeArray
                 = nullptr,
                 const TrackingVolumeVector cVolumeVector = {},
                 const DetachedVolumeVector dVolumeVector = {},
                 const std::string&         volumeName    = "undefined");

  /// Copy Constructor with a shift
  ///
  /// @param tvol is the source tracking volume
  /// @param shift is the additional shift applied after copying
  /// @param volumeName is a string identifier
  ///
  /// @return shared pointer to a new TrackingVolume
  TrackingVolume(const TrackingVolume& tvol,
                 const Transform3D&    shift,
                 const std::string&    volumeName = "undefined");

private:
  /// Create Boundary Surface
  void
  createBoundarySurfaces();

  /// method to synchronize the layers with potentially updated volume bounds:
  /// - adapts the layer dimensions to the new volumebounds + envelope
  ///
  /// @param envelope is the clearance between volume boundary and layer
  void
  synchronizeLayers(double envelope = 1.) const;

  /// close the Geometry, i.e. set the TDD_ID
  ///
  /// @param volumeID is the geometry id of the volume
  ///        as calculated by the TrackingGeometry
  /// @param volumeMap is a map to find the a volume
  ///        by a given name
  void
  closeGeometry(GeometryID& volumeID,
                std::map<std::string, const TrackingVolume*>& volumeMap) const;

  /// interlink the layers in this TrackingVolume
  void
  interlinkLayers();

  /// Forbidden copy constructor - deleted
  TrackingVolume(const TrackingVolume&) = delete;

  /// Forbidden assignment - deleted
  TrackingVolume&
  operator=(const TrackingVolume&)
      = delete;

  /// The Material the TrackingVolume consists of
  std::shared_ptr<Material> m_material;

  /// Remember the mother volume
  mutable const TrackingVolume* m_motherVolume;

  // the boundary surfaces
  mutable std::vector<TrackingVolumeBoundaryPtr> m_boundarySurfaces;

  ///(a) static configuration ordered by Binned arrays
  /// static layers
  mutable std::unique_ptr<const LayerArray> m_confinedLayers;

  /// Array of Volumes inside the Volume
  mutable std::shared_ptr<const TrackingVolumeArray> m_confinedVolumes;

  /// (b)  non-static setups
  /// detacathd
  const DetachedVolumeVector m_confinedDetachedVolumes;
  /// confined dense
  const TrackingVolumeVector m_confinedDenseVolumes;
  /// confined arbitrary
  const LayerVector m_confinedArbitraryLayers;

  /// Volumes to glue Volumes from the outside
  mutable GlueVolumesDescriptor* m_glueVolumeDescriptor;

  /// The Signature done by the GeometryBuilder
  mutable GeometrySignature m_geometrySignature;

  /// The gometry type for the navigation schema
  mutable GeometryType m_geometryType;

  //// Volume name for debug reasons & screen output
  std::string m_name;

  /// color code for displaying
  mutable unsigned int m_colorCode;
};

inline const std::string&
TrackingVolume::volumeName() const
{
  return m_name;
}

inline const LayerArray*
TrackingVolume::confinedLayers() const
{
  return m_confinedLayers.get();
}

inline const LayerVector
TrackingVolume::confinedArbitraryLayers() const
{
  return m_confinedArbitraryLayers;
}

inline std::shared_ptr<const TrackingVolumeArray>
TrackingVolume::confinedVolumes() const
{
  return m_confinedVolumes;
}

inline std::shared_ptr<const TrackingVolumeArray>
TrackingVolume::confinedVolumesSharedPtr() const
{
  return m_confinedVolumes;
}

inline const DetachedVolumeVector
TrackingVolume::confinedDetachedVolumes() const
{
  return m_confinedDetachedVolumes;
}

inline const TrackingVolumeVector
TrackingVolume::confinedDenseVolumes() const
{
  return m_confinedDenseVolumes;
}

template <class T>
bool
TrackingVolume::onVolumeBoundary(const T& pars) const
{
  // get the associated Surface
  const Surface* pSurface  = &pars.associatedSurface();
  auto&          bSurfaces = boundarySurfaces();
  // fast loop pointer comparison of the surfaces
  for (auto& bsIter : bSurfaces) {
    const BoundarySurfaceT<TrackingVolume>* bSurface = bsIter.get();
    // pointer of the parameter surface is identical with one of the boundary
    // surface pointers
    if (pSurface == &bSurface->surfaceRepresentation()) return true;
  }
  // slow loop - checking the onSurface (does pointer comparison as well)
  for (auto& bsIter : bSurfaces) {
    const BoundarySurfaceT<TrackingVolume>* bSurface = bsIter.get();
    // pointer of the parameter surface is identical with one of the boundary
    // surface pointers
    if (bSurface->onBoundary(pars)) return true;
  }
  // could not find an onSurface
  return false;
}

template <class T>
std::vector<LayerIntersection<T>>
TrackingVolume::layerCandidatesOrdered(const Layer*         sLayer,
                                       const Layer*         eLayer,
                                       const T&             pars,
                                       PropDirection        pDir,
                                       const BoundaryCheck& bcheck,
                                       bool                 resolveMaterial,
                                       bool resolveSensitive) const
{
  // get position and momentum from the parameters
  const Vector3D& gp = pars.position();
  const Vector3D& gm = pars.momentum();

  // the layer intersections
  std::vector<LayerIntersection<T>> lIntersections;
  // assign the direction
  const Vector3D& dir
      = (pDir == alongMomentum ? gm.unit() : Vector3D(-1 * gm.unit()));
  // the confinedLayers
  if (m_confinedLayers) {
    // cache the longest path length to avoid punch-through to the other side
    Intersection   sLayerIntersection(Vector3D(0., 0., 0), 0., true, 0.);
    const Surface* sLayerSurface   = 0;
    double         validPathLength = 0.;
    // - get compatible layers back from the LayerArray simply because of the
    // binning
    // start layer given or not - test layer
    const Layer* tLayer = sLayer ? sLayer : associatedLayer(gp);
    if (tLayer) {
      do {
        // collect material or sensitive layers, usually ignore navigation
        // layers, but
        // always provide the final layer for the navigation stop
        if ((resolveMaterial && tLayer->hasMaterial())
            || (resolveSensitive && tLayer->hasSensitive())
            || tLayer == eLayer) {
          // taking the layer representation, not the approach surface, latter
          // is resolved in the navigation stream
          const Surface& tSurface = tLayer->surfaceRepresentation();
          // intersect the layer @todo should probably be surface on approach
          Intersection lIntersection
              = tSurface.intersectionEstimate(gp, dir, true, bcheck);
          // (a) if the current layer is NOT the start layer - intersection is
          // ok
          if (tLayer != sLayer && lIntersection.valid) {
            lIntersections.push_back(LayerIntersection<T>(
                lIntersection, tLayer, &tSurface, 0, pDir));
            validPathLength = lIntersection.pathLength;
          } else if (tLayer == sLayer) {
            // (b) the current layer is the start layer - we need to cache it
            // and check with the path length
            //     this avoids potential punch-through to other side of
            sLayerIntersection = lIntersection;
            sLayerSurface      = &tSurface;
          } else if (tLayer == eLayer) {
            // (c) it is the end layer after all - provide it and break the loop
            lIntersections.push_back(LayerIntersection<T>(
                lIntersection, tLayer, &tSurface, 0, pDir));
            break;
          }
        }
        // move to next one or break because you reached the end layer
        tLayer = (tLayer == eLayer) ? nullptr : tLayer->nextLayer(gp, dir);
      } while (tLayer);
    }

    // final check for compatibility of the start layer in order to avoid
    // punch-through
    if (sLayer && sLayerIntersection.valid
        && sLayerIntersection.pathLength < validPathLength)
      lIntersections.push_back(LayerIntersection<T>(
          sLayerIntersection, sLayer, sLayerSurface, 0, pDir));
  }
  // and the arbitraray layers
  if (!m_confinedArbitraryLayers.empty()) {
    // loop over the layers and intersect them
    for (auto& layer : m_confinedArbitraryLayers) {
      // intersections
      Intersection lIntersection
          = layer->surfaceRepresentation().intersectionEstimate(
              gp, dir, true, bcheck);
      if (lIntersection.valid)
        lIntersections.push_back(
            LayerIntersection<T>(lIntersection,
                                 layer.get(),
                                 &(layer->surfaceRepresentation()),
                                 0,
                                 pDir));
    }
  }

  // sort them accordingly to the path length
  std::sort(lIntersections.begin(), lIntersections.end());
  // and return
  return lIntersections;
}

// Returns the boundary surfaces ordered in probability to hit them based on
// straight line intersection @todo change hard-coded default 
template <class T>
std::vector<BoundaryIntersection<T>>
TrackingVolume::boundarySurfacesOrdered(const T& pars, PropDirection pDir) const
{
  // assign the direction
  const Vector3D dir
      = (pDir == alongMomentum ? pars.momentum().unit()
                               : Vector3D(-1 * pars.momentum().unit()));
  // loop over boundarySurfaces and calculate the intersection
  std::vector<BoundaryIntersection<T>> bIntersections;
  auto&                                bSurfaces = boundarySurfaces();
  for (auto& bsIter : bSurfaces) {
    const BoundarySurfaceT<TrackingVolume>* bSurface = bsIter.get();
    Intersection                            bsIntersection
        = bSurface->surfaceRepresentation().intersectionEstimate(
            pars.position(), dir, true, false);
    if (bsIntersection.valid)
      bIntersections.push_back(
          BoundaryIntersection<T>(bsIntersection,
                                  bSurface,
                                  &(bSurface->surfaceRepresentation()),
                                  0,
                                  pDir));
  }
  // and now sort to get the closest
  std::sort(bIntersections.begin(), bIntersections.end());
  // and return
  return bIntersections;
}

inline GeometrySignature
TrackingVolume::geometrySignature() const
{
  return m_geometrySignature;
}

inline GeometryType
TrackingVolume::geometryType() const
{
  return m_geometryType;
}

inline void
TrackingVolume::registerColorCode(unsigned int icolor) const
{
  m_colorCode = icolor;
}

inline unsigned int
TrackingVolume::colorCode() const
{
  return m_colorCode;
}

inline const TrackingVolume*
TrackingVolume::motherVolume() const
{
  return m_motherVolume;
}

inline void
TrackingVolume::setMotherVolume(const TrackingVolume* mvol) const
{
  m_motherVolume = mvol;
}

}  // end of namespace

#endif  // ACTS_DETECTOR_TRACKINGVOLUME_H

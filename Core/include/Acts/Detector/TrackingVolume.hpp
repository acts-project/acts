// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingVolume.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <map>
#include <string>
#include "Acts/Layers/Layer.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometrySignature.hpp"
#include "Acts/Volumes/BoundarySurfaceT.hpp"
#include "Acts/Volumes/Volume.hpp"

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
typedef std::shared_ptr<TrackingVolume>               MutableTrackingVolumePtr;
typedef std::shared_ptr<const DetachedTrackingVolume> DetachedTrackingVolumePtr;

// possible contained
typedef BinnedArray<TrackingVolumePtr>         TrackingVolumeArray;
typedef std::vector<TrackingVolumePtr>         TrackingVolumeVector;
typedef BinnedArray<LayerPtr>                  LayerArray;
typedef std::vector<LayerPtr>                  LayerVector;
typedef std::vector<DetachedTrackingVolumePtr> DetachedVolumeVector;

// full intersection with Layer
using LayerIntersection = FullIntersection<Layer, Surface>;

// full intersection with surface
using BoundaryIntersection
    = FullIntersection<BoundarySurfaceT<TrackingVolume>, Surface>;

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
  /// @param htrans is the global 3D transform to position the volume in space
  /// @param volumeBounds is the description of the volume boundaries
  /// @param containedVolumes are the static volumes that fill this volume
  /// @param volumeName is a string identifier
  ///
  /// @return shared pointer to a new TrackingVolume
  static MutableTrackingVolumePtr
  create(std::shared_ptr<const Transform3D>         htrans,
         VolumeBoundsPtr                            volumeBounds,
         std::shared_ptr<const TrackingVolumeArray> containedVolumes = nullptr,
         const std::string&                         volumeName = "undefined")
  {
    return MutableTrackingVolumePtr(
        new TrackingVolume(htrans, volumeBounds, containedVolumes, volumeName));
  }

  /// Factory constructor for Tracking Volumes with content
  /// - can not be a container volume
  ///
  /// @param htrans is the global 3D transform to position the volume in space
  /// @param volumeBounds is the description of the volume boundaries
  /// @param matprop is are materials of the tracking volume
  /// @param cLayerArray is the confined layer array (optional)
  /// @param cLayerVector is the confined arbitrary layer vector
  /// @param cVolumeVector is the confined arbitrary volume vector
  /// @param dVolumeVector is the confined arbeitrary detached volume vector
  /// @param volumeName is a string identifier
  ///
  /// @return shared pointer to a new TrackingVolume
  static MutableTrackingVolumePtr
  create(std::shared_ptr<const Transform3D> htrans,
         VolumeBoundsPtr                    volumeBounds,
         std::shared_ptr<const Material>    matprop,
         std::unique_ptr<const LayerArray>  cLayerArray   = nullptr,
         const LayerVector                  cLayerVector  = {},
         const TrackingVolumeVector         cVolumeVector = {},
         const DetachedVolumeVector         dVolumeVector = {},
         const std::string&                 volumeName    = "undefined")
  {
    return MutableTrackingVolumePtr(new TrackingVolume(htrans,
                                                       volumeBounds,
                                                       matprop,
                                                       std::move(cLayerArray),
                                                       cLayerVector,
                                                       nullptr,
                                                       cVolumeVector,
                                                       dVolumeVector,
                                                       volumeName));
  }

  /// Return the associated Layer to the global position
  ///
  /// @param gpos is the associated global position
  ///
  /// @return plain pointer to layer object
  const Layer*
  associatedLayer(const Vector3D& gpos) const;

  /// @brief Resolves the volume into (compatible) Layers
  ///
  /// This is the method for the propagator/extrapolator
  /// @param_t parameters_t Type of parameters used for the decomposition
  /// @param_t options_t Type of navigation options object for decomposition
  /// @param_t corrector_t Type of (optional) corrector for surface intersection
  ///
  /// @param parameters The templated parameters for searching
  /// @param options The templated navigation options
  /// @param corrfnc is the corrector struct / function
  ///
  /// @return vector of compatible intersections with layers
  template <typename parameters_t,
            typename options_t,
            typename corrector_t = VoidCorrector>
  std::vector<LayerIntersection>
  compatibleLayers(const parameters_t& parameters,
                   const options_t&    options,
                   const corrector_t&  corrfnc = corrector_t()) const;

  /// @brief Returns the boundary surfaces ordered in probability
  ///
  /// @param_t parameters_t Type of parameters used for the decomposition
  /// @param_t options_t Type of navigation options object for decomposition
  /// @param_t corrector_t Type of (optional) corrector for surface intersection
  ///
  /// @param parameters The templated parameters for searching
  /// @param options The templated navigation options
  /// @param corrfnc is the corrector struct / function
  ///
  /// @return is the templated boundary intersection
  template <typename parameters_t,
            typename options_t,
            typename corrector_t = VoidCorrector>
  std::vector<BoundaryIntersection>
  compatibleBoundaries(const parameters_t& parameters,
                       const options_t&    options,
                       const corrector_t&  corrfnc = corrector_t()) const;

  /// LEGACY method for old ExtrapolationEngine -------------------------------
  ///
  /// --------- to be deprecated with release 0.07.00
  /// - it will be replace by compatibleLayers()
  ///
  /// Return the material layers ordered based on straight line intersections:
  ///
  /// - startObject and endLayer are included in the list
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
  template <typename parameters_t>
  std::vector<LayerIntersection>
  layerCandidatesOrdered(const Layer*         sLayer,
                         const Layer*         eLayer,
                         const parameters_t&  parameters,
                         NavigationDirection  pDir             = forward,
                         const BoundaryCheck& bcheck           = true,
                         bool                 resolveSensitive = true,
                         bool                 resolveMaterial  = true,
                         bool                 resolvePassive   = false) const;

  /// LEGACY method for old ExtrapolationEngine -------------------------------
  ///
  /// --------- to be deprecated with release 0.07.00
  /// - it will be replace by compatibleBoundaries()
  ///
  /// Returns the boundary surfaces ordered in probability to hit them based
  /// on straight line intersection
  ///
  /// @param parameters are the templated tracking parameters
  /// @param pDir is the additional direction presciprion
  ///
  /// @return is the templated boundary intersection
  template <typename parameters_t>
  std::vector<BoundaryIntersection>
  boundarySurfacesOrdered(const parameters_t& parameters,
                          NavigationDirection pDir        = forward,
                          bool                skipCurrent = false) const;

  /// Return the associated sub Volume, returns THIS if no subVolume exists
  ///
  /// @param gpos is the global position associated with that search
  ///
  /// @return plain pointer to associated with the position
  const TrackingVolume*
  trackingVolume(const Vector3D& gpos) const;

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

  /// Glue another tracking volume to this one
  ///  - if common face is set the glued volumes are sharing the boundary, down
  /// to the last navigation volume
  ///
  /// @param bsfMine is the boundary face indicater where to glue
  /// @param neighbor is the TrackingVolume to be glued
  /// @param bsfNeighbor is the boudnary surface of the neighbor
  void
  glueTrackingVolume(BoundarySurfaceFace      bsfMine,
                     MutableTrackingVolumePtr neighbor,
                     BoundarySurfaceFace      bsfNeighbor);

  /// Glue another tracking volume to this one
  ///  - if common face is set the glued volumes are sharing the boundary, down
  /// to the last navigation volume
  ///
  ///
  /// @param bsfMine is the boundary face indicater where to glue
  /// @param neighbors are the TrackingVolumes to be glued
  /// @param bsfNeighbors are the boudnary surface of the neighbors
  void
  glueTrackingVolumes(BoundarySurfaceFace                  bsfMine,
                      std::shared_ptr<TrackingVolumeArray> neighbors,
                      BoundarySurfaceFace                  bsfNeighbors);

  /// Provide a new BoundarySurface from the glueing
  ///
  ///
  /// @param bsfMine is the boundary face indicater where to glue
  /// @param bsSurface is the new boudnary surface
  void
  updateBoundarySurface(
      BoundarySurfaceFace                                     bsfMine,
      std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bsSurface);

  /// Register the outside glue volumes -
  /// ordering is in the TrackingVolume Frame:
  ///  - negativeFaceXY
  ///  - (faces YZ, ZY, radial faces)
  ///  - positiveFaceXY
  ///
  /// @param gvd register a new GlueVolumeDescriptor
  /// @todo update to shared/unique ptr
  void
  registerGlueVolumeDescriptor(GlueVolumesDescriptor* gvd);

  /// Register the outside glue volumes -
  /// ordering is in the TrackingVolume Frame:
  ///  - negativeFaceXY
  ///  - (faces YZ, ZY, radial faces)
  ///  - positiveFaceXY
  GlueVolumesDescriptor&
  glueVolumesDescriptor();

  /// Sign the volume - the geometry builder has to do that
  ///
  /// @param signat is the volume signature
  /// @param geotype is the volume navigation type
  void
  sign(GeometrySignature signat, GeometryType geotype = Static);

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
  registerColorCode(unsigned int icolor);

  /// Get the color code
  unsigned int
  colorCode() const;

  /// Return the MotherVolume - if it exists
  const TrackingVolume*
  motherVolume() const;

  /// Set the MotherVolume
  ///
  /// @param mvol is the mother volume
  void
  setMotherVolume(const TrackingVolume* mvol);

  /// @return a map of all contained detector elements and their corresponding
  /// identifier
  const std::map<Identifier, const DetectorElementBase*>&
  detectorElements() const;

protected:
  /// Default constructor
  TrackingVolume();

  /// Constructor for a container Volume
  /// - vacuum filled volume either as a for other tracking volumes
  ///
  /// @param htrans is the global 3D transform to position the volume in space
  /// @param volumeBounds is the description of the volume boundaries
  /// @param containedVolumes are the static volumes that fill this volume
  /// @param volumeName is a string identifier
  TrackingVolume(
      std::shared_ptr<const Transform3D>               htrans,
      VolumeBoundsPtr                                  volumeBounds,
      const std::shared_ptr<const TrackingVolumeArray> containedVolumes
      = nullptr,
      const std::string& volumeName = "undefined");

  /// Constructor for a full equipped Tracking Volume
  ///
  /// @param htrans is the global 3D transform to position the volume in space
  /// @param volumeBounds is the description of the volume boundaries
  /// @param matprop is are materials of the tracking volume
  /// @param cLayerArray is the confined layer array (optional)
  /// @param cLayerVector is the confined arbitrary layer vector
  /// @param cVolumeArray is the confined volume array
  /// @param cVolumeVector is the confined arbitrary volume vector
  /// @param dVolumeVector is the confined arbeitrary detached volume vector
  /// @param volumeName is a string identifier
  TrackingVolume(std::shared_ptr<const Transform3D> htrans,
                 VolumeBoundsPtr                    volumeBounds,
                 std::shared_ptr<const Material>    matprop,
                 std::unique_ptr<const LayerArray>  cLayerArray  = nullptr,
                 const LayerVector                  cLayerVector = {},
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
  /// @param volumeMap is a map to find the a volume
  ///        by a given name
  /// @param vol is the geometry id of the volume
  ///        as calculated by the TrackingGeometry
  ///
  void
  closeGeometry(std::map<std::string, const TrackingVolume*>& volumeMap,
                size_t& vol);

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
  std::shared_ptr<const Material> m_material;

  /// Remember the mother volume
  const TrackingVolume* m_motherVolume;

  // the boundary surfaces
  std::vector<TrackingVolumeBoundaryPtr> m_boundarySurfaces;

  ///(a) static configuration ordered by Binned arrays
  /// static layers
  std::unique_ptr<const LayerArray> m_confinedLayers;

  /// Array of Volumes inside the Volume
  std::shared_ptr<const TrackingVolumeArray> m_confinedVolumes;

  /// (b)  non-static setups
  /// detacathd
  const DetachedVolumeVector m_confinedDetachedVolumes;
  /// confined dense
  const TrackingVolumeVector m_confinedDenseVolumes;
  /// confined arbitrary
  const LayerVector m_confinedArbitraryLayers;

  /// Volumes to glue Volumes from the outside
  GlueVolumesDescriptor* m_glueVolumeDescriptor;

  /// The Signature done by the GeometryBuilder
  GeometrySignature m_geometrySignature;

  /// The gometry type for the navigation schema
  GeometryType m_geometryType;

  //// Volume name for debug reasons & screen output
  std::string m_name;

  /// color code for displaying
  unsigned int m_colorCode;

  /// map which collects all detector elements and connects them with their
  /// Identfier
  std::map<Identifier, const DetectorElementBase*> m_detectorElements;
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
TrackingVolume::registerColorCode(unsigned int icolor)
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
TrackingVolume::setMotherVolume(const TrackingVolume* mvol)
{
  m_motherVolume = mvol;
}

inline const std::map<Identifier, const DetectorElementBase*>&
TrackingVolume::detectorElements() const
{
  return m_detectorElements;
}

#include "detail/TrackingVolume.ipp"

}  // namespace

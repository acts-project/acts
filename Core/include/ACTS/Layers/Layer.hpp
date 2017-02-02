// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Layer.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_LAYER_H
#define ACTS_DETECTOR_LAYER_H 1

// Core module
#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Utilities/ApproachDescriptor.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/GeometryObject.hpp"
#include "ACTS/Utilities/Intersection.hpp"
#include "ACTS/Volumes/AbstractVolume.hpp"

namespace Acts {

class Surface;
class SurfaceMaterial;
class MaterialProperties;
class BinUtility;
class Volume;
class VolumeBounds;
class TrackingVolume;
class DetachedTrackingVolume;
class ApproachDescriptor;
class ICompatibilityEstimator;

typedef ObjectIntersection<Surface> SurfaceIntersection;

// master typedef
class Layer;
typedef std::shared_ptr<const Layer> LayerPtr;
typedef std::shared_ptr<Layer>       MutableLayerPtr;
typedef std::pair<const Layer*, const Layer*> NextLayers;

///
/// @enum LayerType
///
/// For readability
///
enum LayerType { navigation = -1, passive = 0, active = 1 };

/// @class Layer
///
/// Base Class for a Detector Layer in the Tracking realm.
/// An actual implemented Detector Layer inheriting from this base
/// class has to inherit from a specific type of Surface as well.
/// In addition, a Layer can carry:
///
/// - a SurfaceArray of Surfaces holding the actual detector elements or
/// subSurfaces.
/// - SurfaceMaterial for Surface-based materialUpdates
/// - a pointer to the TrackingVolume (can only be set by such)
/// - an active/passive code :
/// 0      - activ
/// 1      - passive
/// [....] - other
///
/// The search type for compatible surfaces on a layer is
///   [ the higher the number, the faster ]:
/// --------- Layer internal ------------------------------------------------
/// -1     - untested: provide all layer surfaces to the extrapolation engine
///               - does not work with endSurface,
///                 will be increased to 0 if endSurface is given
///               - debug mode only !
///  0     - test all on intersection and provide to the extrapolation engine
///  1     - provide bin surface and registered neighbours and bin mates
///               - does not work with endSurface,
///                 will be increased to 2 if endSurface is given
///  2     - as 1 but with intersection test @todo compatibility test
/// @todo check for update
///  3     - provide bin surface and next bin surfaces (if differ)
///               - does not work with endSurface
///                 will be increased to 4 if endSurface is given
///  4     - as 3 but with intersection test @todo compatibility test
///  5     - whatever the overlap descriptor returns with this
///
///
class Layer : public virtual GeometryObject
{
  /// Declare the TrackingVolume as a friend, to be able to register previous,
  /// next and set the enclosing TrackingVolume
  friend class TrackingVolume;

  /// Declare the DetachedTrackingVolume as a friend to be able to register it
  friend class DetachedTrackingVolume;

public:
  /// Clone at a with a shift - this is the only clone allowed
  ///
  /// @param shift is the additional transform applied after cloning
  virtual LayerPtr
  cloneWithShift(const Transform3D& shift) const = 0;

  /// Destructor
  virtual ~Layer();

  /// Assignment operator - forbidden, layer assignment must not be ambiguous
  ///
  /// @param lay is the source layer for assignment
  Layer&
  operator=(const Layer& lay)
      = delete;

  /// Return the entire SurfaceArray, returns a nullptr if no SurfaceArray
  const SurfaceArray*
  surfaceArray() const;

  /// Transforms the layer into a Surface representation for extrapolation
  /// @note the layer can be hosting many surfaces, but this is the global
  /// one to which one can extrapolate
  virtual const Surface&
  surfaceRepresentation() const = 0;
  
  // Non-const version
  virtual Surface&
  surfaceRepresentation() = 0;

  /// Return the Thickness of the Layer
  /// this is by definition along the normal vector of the surfaceRepresentation
  double
  thickness() const;

  /// templated onLayer() method
  ///
  /// @param parameters are the templated (charged/neutral) on layer check
  /// @param bcheck is the boundary check directive
  ///
  /// @return boolean that indicates success of the operation
  template <class T>
  bool
  onLayer(const T& parameters, const BoundaryCheck& bcheck = true) const;

  /// geometrical isOnLayer() method
  ///
  /// @note using isOnSurface() with Layer specific tolerance
  /// @param pos is the gobal position to be checked
  /// @param bcheck is the boundary check directive
  ///
  /// @return boolean that indicates success of the operation
  virtual bool
  isOnLayer(const Vector3D& pos, const BoundaryCheck& bcheck = true) const;

  /// Return method for the approach descriptor, can be nullptr
  const ApproachDescriptor*
  approachDescriptor() const;
  
  /// Non-const version
  ApproachDescriptor*
  approachDescriptor();

  ///  Surface seen on approach
  /// for surcfaces without sub structure, this is the surfaceRepresentation
  /// @param gpos is the global position to start the approach from
  /// @param dir is the direction from where to attempt the approach
  /// @param pdir is the direction prescription
  /// @param bcheck is the boundary check directive
  /// @param resolveSubSurfaces is a boolean directive whether to resolve
  /// structure or not
  /// @note reasons for resolving are: collect & find material, collect & find
  /// sensitivex
  /// @param ice is a (future) compatibility estimator that could be used to
  /// modify
  /// the straight line approach
  virtual const SurfaceIntersection
  surfaceOnApproach(const Vector3D&                gpos,
                    const Vector3D&                dir,
                    PropDirection                  pdir,
                    const BoundaryCheck&           bcheck,
                    bool                           resolveSubSurfaces = false,
                    const ICompatibilityEstimator* ice = nullptr) const;

  ///  get compatible surfaces starting from charged parameters
  ///  returns back the compatible surfaces either with straight line
  ///  estimation,
  ///  or (@todo later) with a compatiblityEstimator.
  ///  - if start/end surface are given, surfaces are provided in between (start
  /// & end excluded)
  ///  - the boolean indicates if the surfaces are direction ordered
  ///
  /// @param cSurfaces are the retrun surface intersections
  /// @param pars are the (charged) track parameters for the search
  /// @param pdir is the propagation direction prescription
  /// @param bcheck is the boundary check directive
  /// @param collectSensitive is the prescription to find the sensitive surfaces
  /// @param collectPassive is the prescription to find all passive surfaces
  /// @param searchType is the level of depth for the search
  /// @param startSurface is an (optional) start surface for the search:
  /// excluded in return
  /// @param endSurface is an (optional) end surface for the search: excluded in
  /// return
  /// @param ice is a (future) compatibility estimator that could be used to
  /// modify
  /// the straight line approach
  ///
  /// @return boolean that indicates if a compatible surface exists at all
  virtual bool
  compatibleSurfaces(std::vector<SurfaceIntersection>& cSurfaces,
                     const TrackParameters&            pars,
                     PropDirection                     pdir,
                     const BoundaryCheck&              bcheck,
                     bool                              collectSensitive,
                     bool                              collectPassive,
                     int                               searchType,
                     const Surface*                    startSurface = nullptr,
                     const Surface*                    endSurface   = nullptr,
                     const ICompatibilityEstimator*    ice = nullptr) const;

  ///  get compatible surfaces starting from neutral parameters
  ///  returns back the compatible surfaces either with straight line
  ///  estimation,
  ///  or (@todo later) with a compatiblityEstimator.
  ///  - if start/end surface are given, surfaces are provided in between (start
  /// & end excluded)
  ///  - the boolean indicates if the surfaces are direction ordered
  ///
  /// @param cSurfaces are the retrun surface intersections
  /// @param pars are the (charged) track parameters for the search
  /// @param pdir is the propagation direction prescription
  /// @param bcheck is the boundary check directive
  /// @param collectSensitive is the prescription to find the sensitive surfaces
  /// @param collectPassive is the prescription to find all passive surfaces
  /// @param searchType is the level of depth for the search
  /// @param startSurface is an (optional) start surface for the search:
  /// excluded in return
  /// @param endSurface is an (optional) end surface for the search: excluded in
  /// return
  /// @param ice is a (future) compatibility estimator that could be used to
  /// modify
  /// the straight line approach
  ///
  /// @return boolean that indicates if a compatible surface exists at all
  virtual bool
  compatibleSurfaces(std::vector<SurfaceIntersection>& cSurfaces,
                     const NeutralParameters&          pars,
                     PropDirection                     pdir,
                     const BoundaryCheck&              bcheck,
                     bool                              collectSensitive,
                     bool                              collectPassive,
                     int                               searchType,
                     const Surface*                    startSurface = nullptr,
                     const Surface*                    endSurface   = nullptr,
                     const ICompatibilityEstimator*    ice = nullptr) const;

  /// Has sub-structure method:
  /// @note sub-structure depending on :
  ///   (a) only when required to resolve sub surfaces for sensitive hits
  ///   (b) also material is ordered with sub structure
  ///
  /// @param resolveSensitive is a directive whether
  /// one should look for sensitive surfaces in the surface array
  ///
  /// @return bollean that indicates if sub structure exitst
  virtual bool
  hasSubStructure(bool resolveSensitive = false) const;

  //  Boolean check method if layer has material:
  // - checks if any of the layer surfaces has material:
  // - can be approach surfaces or layer surface
  ///
  /// @return bollean that indicates if material exists
  virtual bool
  hasMaterial() const;

  /// method returning a pointer to the surface material of the layer
  /// @note not the layer itself holds the material but one of its
  /// approachsurfaces
  ///
  /// @return the surface material
  const SurfaceMaterial*
  material() const;

  /// Returns a pointer to the surface which carries the surface material
  /// @note can be either inner, outer BoundarySurface or the surface
  /// representation (central)
  /// return nullptr if layer has not material
  ///
  /// @todo remove this concept
  const Surface*
  materialSurface() const;
  
  // Non-const version
  Surface*
  materialSurface();

  /// Boolean check method if layer has sensitive surfaces
  /// @note checks if a surfaceArray is present
  ///
  /// @todo reemove this concept
  virtual bool
  hasSensitive() const;

  /// Fast navigation to next layer
  ///
  /// @param pos is the start position for the search
  /// @param dir is the direction for the search
  ///
  /// @return the pointer to the next layer
  const Layer*
  nextLayer(const Vector3D& pos, const Vector3D& dir) const;

  /// get the confining TrackingVolume
  ///
  /// @return the pointer to the enclosing volume
  const TrackingVolume*
  enclosingTrackingVolume() const;

  /// get the confining DetachedTrackingVolume
  ///
  /// @return the pointer to the detached volume
  const DetachedTrackingVolume*
  enclosingDetachedTrackingVolume() const;

  /// register Volume associated to the layer - if you want to do that by hand
  ///
  /// @param avol is the provided volume
  void
  registerRepresentingVolume(const AbstractVolume* avol);

  ///  return the abstract volume that represents the layer
  ///
  /// @return the representing volume of the layer
  const AbstractVolume*
  representingVolume() const;

  /// return the LayerType
  LayerType
  layerType() const;

protected:
  /// Default Constructor
  Layer();

  /// Copy Constructor
  ///
  /// @param lay is the source layer for the copye
  Layer(const Layer& lay);

  /// Constructor with pointer to SurfaceArray (passing ownership)
  ///
  /// @param surfaceArray is the array of sensitive surfaces
  /// @param thickness is the normal thickness of the Layer
  /// @param ad oapproach descriptor
  /// @param ltype is the layer type if active or passive
  Layer(std::unique_ptr<SurfaceArray>       surfaceArray,
        double                              thickness = 0.,
        std::unique_ptr<ApproachDescriptor> ad        = nullptr,
        LayerType                           ltype     = passive);

  /// get compatible surfaces starting from charged parameters - forward call
  /// from explicit methods
  ///
  /// @param cSurfaces are the retrun surface intersections
  /// @param pars are the (charged) track parameters for the search
  /// @param pdir is the propagation direction prescription
  /// @param bcheck is the boundary check directive
  /// @param collectSensitive is the prescription to find the sensitive surfaces
  /// @param collectPassive is the prescription to find all passive surfaces
  /// @param searchType is the level of depth for the search
  /// @param startSurface is an (optional) start surface for the search:
  /// excluded in return
  /// @param endSurface is an (optional) end surface for the search: excluded in
  /// return
  /// @param ice is a (future) compatibility estimator that could be used to
  /// modify
  /// the straight line approach
  ///
  /// @return boolean that indicates if a compatible surface exists at all
  template <class T>
  bool
  getCompatibleSurfaces(std::vector<SurfaceIntersection>& cSurfaces,
                        const T&                          pars,
                        PropDirection                     pdir,
                        const BoundaryCheck&              bcheck,
                        bool                              collectSensitive,
                        bool                              collectPassive,
                        int                               searchType,
                        const Surface*                 startSurface = nullptr,
                        const Surface*                 endSurface   = nullptr,
                        const ICompatibilityEstimator* ice = nullptr) const;

  /// test compatible surface - checking directly for intersection & collection
  ///
  /// geometrical test compatible surface method
  /// @param cSurfaces are the retrun surface intersections
  /// @param surface is the parameter surface
  /// @todo how is this different from the start surface
  /// @param gpos is the resolved global position
  /// @param dir is themomentum direction
  /// @param pdir is the propagation direction prescription
  /// @param bcheck is the boundary check directive
  /// @param maxPathLength is the maximal path length allowed to the surface
  /// @param collectSensitive is the prescription to find the sensitive surfaces
  /// @param collectPassive is the prescription to find all passive surfaces
  /// @param intersectionTest is a boolean idicating if intersection is done
  /// @param startSurface is an (optional) start surface for the search:
  /// excluded in return
  /// @param endSurface is an (optional) end surface for the search: excluded in
  /// return
  /// @param ice is a (future) compatibility estimator that could be used to
  /// modify
  /// the straight line approach
  ///
  /// @return boolean that indicates if a compatible surface exists at all
  void
  testCompatibleSurface(std::vector<SurfaceIntersection>& cSurfaces,
                        const Surface&                    surface,
                        const Vector3D&                   gpos,
                        const Vector3D&                   dir,
                        PropDirection                     pdir,
                        const BoundaryCheck&              bcheck,
                        double                            maxPathLength,
                        bool                              collectSensitive,
                        bool                              collectPassive,
                        bool                              intersectionTest,
                        const Surface*                 startSurface = nullptr,
                        const Surface*                 endSurface   = nullptr,
                        const ICompatibilityEstimator* ice = nullptr) const;

  ///  private method to set enclosing TrackingVolume, called by friend class
  /// only
  ///  optionally, the layer can be resized to the dimensions of the
  /// TrackingVolume
  ///  - Bounds of the Surface are resized
  ///  - MaterialProperties dimensions are resized
  ///  - SubSurface array boundaries are NOT resized
  ///
  /// @param tvol is the tracking volume the layer is confined
  void
  encloseTrackingVolume(const TrackingVolume& tvol);

  ///  private method to set the enclosed detached TV,
  /// called by friend class only
  ///
  /// @param dtvol is the detached tracking volume the layer is confined
  void
  encloseDetachedTrackingVolume(const DetachedTrackingVolume& dtvol);

  /// the previous Layer according to BinGenUtils
  NextLayers m_nextLayers;
  /// A binutility to find the next layer
  const BinUtility* m_nextLayerUtility;  // @TODO check if this is needed

  /// SurfaceArray on this layer Surface
  std::unique_ptr<SurfaceArray> m_surfaceArray;
  /// thickness of the Layer
  double m_layerThickness;
  /// descriptor for surface on approach
  std::unique_ptr<ApproachDescriptor> m_approachDescriptor;
  /// the enclosing TrackingVolume
  const TrackingVolume* m_enclosingTrackingVolume;
  /// the eventual enclosing detached Tracking volume
  const DetachedTrackingVolume* m_enclosingDetachedTrackingVolume;
  /// Representing Volume
  /// can be used as appraoch suraces
  const AbstractVolume* m_representingVolume;
  /// make a passive/active divisio
  LayerType m_layerType;
  /// pointer to the approach surface carrying the material
  /// nullptr if no material set
  /// @todo remove this concept
  const Surface* m_materialSurface;

private:
  /// Private helper method to close the geometry
  /// @param layerID is the geometry id of the volume
  ///                as calculated by the TrackingGeometry
  void
  closeGeometry(const GeometryID& layerID);
};

inline const SurfaceArray*
Layer::surfaceArray() const
{
  return m_surfaceArray.get();
}

inline double
Layer::thickness() const
{
  return m_layerThickness;
}

inline LayerType
Layer::layerType() const
{
  return m_layerType;
}

inline const TrackingVolume*
Layer::enclosingTrackingVolume() const
{
  return m_enclosingTrackingVolume;
}

inline void
Layer::encloseTrackingVolume(const TrackingVolume& tvol)
{
  m_enclosingTrackingVolume = &(tvol);
}

inline const DetachedTrackingVolume*
Layer::enclosingDetachedTrackingVolume() const
{
  return m_enclosingDetachedTrackingVolume;
}

inline void
Layer::encloseDetachedTrackingVolume(const DetachedTrackingVolume& tvol)
{
  m_enclosingDetachedTrackingVolume = &(tvol);
}

inline const AbstractVolume*
Layer::representingVolume() const
{
  return m_representingVolume;
}

inline const Layer*
Layer::nextLayer(const Vector3D& gp, const Vector3D& mom) const
{
  // no binutility -> no chance to find out the direction
  if (!m_nextLayerUtility) return nullptr;
  return (m_nextLayerUtility->nextDirection(gp, mom) < 0) ? m_nextLayers.first
                                                          : m_nextLayers.second;
}

inline void
Layer::registerRepresentingVolume(const AbstractVolume* theVol)
{
  delete m_representingVolume;
  m_representingVolume = theVol;
}

/// Layers are constructedd with shared_ptr factories, hence the layer array is
/// describes as:
typedef BinnedArray<LayerPtr> LayerArray;

}  // end of namespace

#include "ACTS/Layers/detail/Layer.ipp"

#endif  // ACTS_DETECTOR_LAYER_H

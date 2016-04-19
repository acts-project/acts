///////////////////////////////////////////////////////////////////
// TrackingVolume.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_TRACKINGVOLUME_H
#define ACTS_DETECTOR_TRACKINGVOLUME_H 1

// Core module
#include "ACTS/Utilities/AlgebraDefinitions.h"
// Geometry module
#include "ACTS/Surface.h"
#include "ACTS/Surfaces/BoundaryCheck.h"
#include "ACTS/Volumes/Volume.h"
//#include "ACTS/Volumes/AbstractVolume.h"
#include "ACTS/Volumes/BoundarySurface.h"
#include "ACTS/Layer.h"
#include "ACTS/Material/Material.h"
#include "ACTS/GeometryUtils/BinnedArray.h"
#include "ACTS/GeometryUtils/GeometrySignature.h"
// EventData module
#include "ACTS/Utilities/PropDirection.h"
//STL
#include <string>

namespace Acts {

  // classes
  class TrackingVolume;
  class DetachedTrackingVolume;
  class GlueVolumesDescriptor;
  class VolumeBounds;

  // master typedefs
  typedef std::shared_ptr<const TrackingVolume>          TrackingVolumePtr;
  typedef std::shared_ptr< const DetachedTrackingVolume> DetachedTrackingVolumePtr;

  // possible contained
  typedef BinnedArray< TrackingVolumePtr >         TrackingVolumeArray;
  typedef std::vector< TrackingVolumePtr >         TrackingVolumeVector;
  typedef BinnedArray< LayerPtr >                  LayerArray;
  typedef std::vector< LayerPtr >                  LayerVector;
  typedef std::vector< DetachedTrackingVolumePtr > DetachedVolumeVector;

  template <class T> using LayerIntersection     = FullIntersection< Layer, Surface, T >;
  template <class T> using BoundaryIntersection  = FullIntersection< BoundarySurface<TrackingVolume>, Surface, T >;

  /**
   @class TrackingVolume

   Full Volume description used in Tracking,
   it inherits from Volume to get the geometrical structure.

       A TrackingVolume at navigation level can provide the (layer) material information / internal navigation with in
       5 different ways:

           --- a) Static confinement of Layers
           --- b) detached sub volumes
           --- b) unordered (arbitrarily oriented) layers
           --- d) unordered sub volumes
           --- e) unordered layers AND unordered subvolumes

      The TrackingVolume can also be a simple container of other TrackingVolumes

   In addition it is capable of holding a subarray of Layers and TrackingVolumes.

   @author Andreas.Salzburger@cern.ch
   */

   class TrackingVolume : public Volume,
                          public Material {

    public:
      /** Destructor */
      ~TrackingVolume();

      /** Factory constructor for a conatiner TrackingVolume - by definition a Vacuum volume */
      static TrackingVolumePtr create(std::shared_ptr<Transform3D> htrans,
                                      VolumeBoundsPtr volumeBounds,
                                      std::shared_ptr<const TrackingVolumeArray> containedVolumes = nullptr,
                                      const std::string& volumeName="undefined")
        { return TrackingVolumePtr(new TrackingVolume(htrans, volumeBounds, containedVolumes, volumeName)); }

      /** Factory constructor for Tracking Volumes with content - can not be a container volume */
      static TrackingVolumePtr create(std::shared_ptr<Transform3D> htrans,
                                      VolumeBoundsPtr volumeBounds,
                                      const Material& matprop,
                                      const LayerArray* cLayerArray = nullptr,
                                      const LayerVector* cLayerVector = nullptr,
                                      const TrackingVolumeVector* cVolumeVector = nullptr,
                                      const DetachedVolumeVector* dVolumeVector = nullptr,
                                      const std::string& volumeName="undefined")
        { return TrackingVolumePtr(new TrackingVolume(htrans,
                                                      volumeBounds,
                                                      matprop,
                                                      cLayerArray,
                                                      cLayerVector,
                                                      nullptr,
                                                      cVolumeVector,
                                                      dVolumeVector,
                                                      volumeName)); }
      /** Factory constructor with a shift */
      static TrackingVolumePtr create(const TrackingVolume& trVol,
                                      const Transform3D& shift,
                                      const std::string& volumeName="undefined" )
       { return TrackingVolumePtr(new TrackingVolume(trVol,shift,volumeName)); }

      /** Virtual constructor clone at new position */
      TrackingVolumePtr cloneWithShift(const Transform3D& shift, const std::string& volumeName="undefined") const
        { return TrackingVolume::create(*this,shift,volumeName); }

      /** Return the associated Layer */
      const Layer* associatedLayer(const Vector3D& gp) const;

      /** Return the material layers ordered based on straight line intersections:
          - startLayer and endLayer are included in the list */
      template <class T> std::vector< LayerIntersection<T> >
                          layerCandidatesOrdered(const Layer* sLayer,
                                                 const Layer* eLayer,
                                                 const T& parameters,
                                                 PropDirection pDir        = alongMomentum,
                                                 const BoundaryCheck& bchk = true,
                                                 bool resolveMaterial      = true,
                                                 bool resolveSubSurfaces   = false) const;

      /** Return the associated sub Volume, returns THIS if no subVolume exists */
      const TrackingVolume* associatedSubVolume(const Vector3D& gp) const;

      /** Return the next volume along the navigation stream */
      const TrackingVolume* nextVolume(const Vector3D& gp, const Vector3D& dir, PropDirection pDir = alongMomentum) const;

      /** Return the dynamically created vector of detached sub volumes */
      const DetachedVolumeVector* assocDetachedSubVolumes(const Vector3D& gp, double tol) const;

      /** Return the confined static layer array - if it exists */
      const LayerArray* confinedLayers() const;

      /** Return the arbitrary layer array - if it exists */
      const LayerVector* confinedArbitraryLayers() const;

      /** Return the confined volumes of this container array - if it exists */
      const TrackingVolumeArray* confinedVolumes() const;

      /** Return the confind volume array as a shared pointer - for glueing */
      std::shared_ptr<const TrackingVolumeArray> confinedVolumesSharedPtr() const;

      /** Return detached subVolumes - not the ownership */
      const DetachedVolumeVector* confinedDetachedVolumes() const;

      /** Return unordered subVolumes - not the ownership */
      const TrackingVolumeVector* confinedDenseVolumes() const;

      /** Returns the VolumeName - for debug reason, might be depreciated later */
      const std::string& volumeName() const;

      /** Method to return the BoundarySurfaces */
      const std::vector< std::shared_ptr<const BoundarySurface<TrackingVolume> > >& boundarySurfaces() const;

      /** Returns the boundary surfaces ordered in probability to hit them based on straight line intersection */
      template <class T> std::vector< BoundaryIntersection<T> >
                         boundarySurfacesOrdered(const T& parameters,
                                                 PropDirection pDir = alongMomentum,
                                                 bool startOffBoundary = false) const;

      /** show if you are on a boundary surface */
      template <class T> bool onVolumeBoundary(const T& pars ) const;


      /** glue another tracking volume to this one
          - if common face is set the glued volumes are sharing the boundary, down to the last navigation volume
      */
      void glueTrackingVolume(BoundarySurfaceFace bsfMine,
                              TrackingVolumePtr neighbor,
                              BoundarySurfaceFace bsfNeighbor) const;

      /** glue another tracking volume to this one
          - if common face is set the glued volumes are sharing the boundary, down to the last navigation volume
      */
      void glueTrackingVolumes(BoundarySurfaceFace bsfMine,
                               std::shared_ptr<const TrackingVolumeArray> neighbors,
                               BoundarySurfaceFace bsfNeighbors) const;

      /** provide a new BoundarySurface from the glueing */
      void updateBoundarySurface(BoundarySurfaceFace bf, std::shared_ptr<const BoundarySurface<TrackingVolume> > bs) const;

      /** Register the outside glue volumes -
          ordering is in the TrackingVolume Frame:
           - negativeFaceXY
           - (faces YZ, ZY, radial faces)
           - positiveFaceXY */

      void registerGlueVolumeDescriptor(GlueVolumesDescriptor* gvd) const;

      /** Register the outside glue volumes -
          ordering is in the TrackingVolume Frame:
           - negativeFaceXY
           - (faces YZ, ZY, radial faces)
           - positiveFaceXY */
      const GlueVolumesDescriptor& glueVolumesDescriptor() const;

      /** sign the volume - the geometry builder has to do that */
      void sign(GeometrySignature signat, GeometryType gtype = Static) const;

      /** return the Signature */
      GeometrySignature geometrySignature() const;

      /** return the Signature */
      GeometryType geometryType() const;

      /** Register the color code */
      void registerColorCode(unsigned int icolor) const;

      /** Get the color code */
      unsigned int colorCode() const;

      /** Return the MotherVolume - if it exists */
      const TrackingVolume* motherVolume() const;

      /** Return the MotherVolume - if it exists */
      void setMotherVolume(const TrackingVolume* mvol) const;

      /** add Material */
      void addMaterial( const Material& mat, float fact=1. ) const;

    protected:
      /** Default constructor */
      TrackingVolume();

      /** Constructor for a container Volume
      - vacuum filled volume either as a for other tracking volumes */
      TrackingVolume(std::shared_ptr<Transform3D> htrans,
                     VolumeBoundsPtr volbounds,
                     const std::shared_ptr<const TrackingVolumeArray> subVolumes = nullptr,
                     const std::string& volumeName="undefined");

      /** Constructor for a full equipped Tracking Volume  */
      TrackingVolume(std::shared_ptr<Transform3D> htrans,
                     VolumeBoundsPtr volbounds,
                     const Material& matprop,
                     const LayerArray* cLayerArray = nullptr,
                     const LayerVector* cLayerVector = nullptr,
                     const TrackingVolumeArray* cVolumeArray = nullptr,
                     const TrackingVolumeVector* cVolumeVector = nullptr,
                     const DetachedVolumeVector* dVolumeVector = nullptr,
                     const std::string& volumeName="undefined");

      /** Copy Constructor with a shift  */
      TrackingVolume(const TrackingVolume& tvol,
                     const Transform3D& shift,
                     const std::string& volumeName="undefined");

    private:
      /** Create Boundary Surface */
      void createBoundarySurfaces();

      /** method to synchronize the layers with potentially updated volume bounds:
          - adapts the layer dimensions to the new volumebounds + envelope */
      void synchronizeLayers(double envelope = 1.) const;

      /** interlink the layers in this TrackingVolume */
      void interlinkLayers();

      /** Forbidden copy constructor */
      TrackingVolume(const TrackingVolume&): Volume(), Material() {}

      /** Forbid assignment */
      TrackingVolume &operator=(const TrackingVolume&) { return *this; }

      mutable const TrackingVolume*                                                    m_motherVolume;            //!< mother volume of this volume

      mutable std::vector< std::shared_ptr<const BoundarySurface<TrackingVolume> > >*  m_boundarySurfaces;        //!< boundary Surfaces

      //(a) static configuration ordered by Binned arrays
      mutable const LayerArray*                                                        m_confinedLayers;          //!< Array of Layers inside the Volume
      mutable std::shared_ptr<const TrackingVolumeArray>                               m_confinedVolumes;         //!< Array of Volumes inside the Volume
      //(b)  non-static setups
      mutable const DetachedVolumeVector*                                              m_confinedDetachedVolumes; //!< Detached subvolumes
      mutable const TrackingVolumeVector*                                              m_confinedDenseVolumes;    //!< Unordered subvolumes
      mutable const LayerVector*                                                       m_confinedArbitraryLayers; //!< Unordered Layers inside the Volume

      mutable GlueVolumesDescriptor*                                                   m_glueVolumeDescriptor;      //!< Volumes to glue Volumes from the outside

      mutable GeometrySignature                                                        m_geometrySignature;       //!< The Signature done by the GeometryBuilder
      mutable GeometryType                                                             m_geometryType;            //!< defines how the extrapolation type

      std::string                                                                      m_name;                    //!< Volume name for debug reasons
      mutable unsigned int                                                             m_colorCode;               //!< Color code for displaying

  };

  inline const std::string& TrackingVolume::volumeName() const { return m_name; }

  inline const LayerArray*  TrackingVolume::confinedLayers() const
  { return m_confinedLayers; }

  inline const LayerVector* TrackingVolume::confinedArbitraryLayers() const
  { return m_confinedArbitraryLayers; }

  inline const TrackingVolumeArray* TrackingVolume::confinedVolumes() const
  { return m_confinedVolumes.get(); }

  inline std::shared_ptr<const TrackingVolumeArray> TrackingVolume::confinedVolumesSharedPtr() const
  { return m_confinedVolumes; }

  inline const DetachedVolumeVector* TrackingVolume::confinedDetachedVolumes() const
  { return m_confinedDetachedVolumes; }

  inline const TrackingVolumeVector* TrackingVolume::confinedDenseVolumes() const
  { return m_confinedDenseVolumes; }

  template <class T> bool TrackingVolume::onVolumeBoundary(const T& pars) const {
      // get the associated Surface
      const Surface* pSurface = &pars.associatedSurface();
      auto& bSurfaces = boundarySurfaces();
      // fast loop pointer comparison of the surfaces
      for (auto& bsIter : bSurfaces ){
          const BoundarySurface<TrackingVolume>* bSurface = bsIter.get();
          // pointer of the parameter surface is identical with one of the boundary surface pointers
          if (pSurface == &bSurface->surfaceRepresentation()) return true;
      }
      // slow loop - checking the onSurface (does pointer comparison as well)
      for (auto& bsIter : bSurfaces ){
          const BoundarySurface<TrackingVolume>* bSurface = bsIter.get();
          // pointer of the parameter surface is identical with one of the boundary surface pointers
          if (bSurface->onBoundary(pars)) return true;
      }
      // could not find an onSurface
      return false;
  }


  /** Return the material layers ordered based on straight line intersections
      - start and end layer are always part of it */
  template <class T> std::vector< LayerIntersection<T> >
                       TrackingVolume::layerCandidatesOrdered(const Layer* sLayer,
                                                              const Layer* eLayer,
                                                              const T& pars,
                                                              PropDirection pDir,
                                                              const BoundaryCheck& bchk,
                                                              bool resolveMaterial,
                                                              bool resolveSensitive) const
  {
      // get position and momentum from the parameters
      const Vector3D& gp = pars.position();
      const Vector3D& gm = pars.momentum();

      // the layer intersections
      std::vector< LayerIntersection<T> > lIntersections;
      // assign the direction
      const Vector3D& dir = ( pDir == alongMomentum ? gm.unit() : Vector3D(-1*gm.unit()));
      // the confinedLayers
      if (m_confinedLayers){

          // cache the longest path length to avoid punch-through to the other side
          Intersection     sLayerIntersection(Vector3D(0.,0.,0),0.,true,0.);
          const Surface*   sLayerSurface  = 0;
          double validPathLength = 0.;
          // - get compatible layers back from the LayerArray simply because of the binning
          // start layer given or not - test layer
          const Layer* tLayer = sLayer ? sLayer : associatedLayer(gp);
          if (tLayer){
              do {
                  // collect material or sensitive layers, usually ignore navigation layers, but
                  // always provide the final layer for the navigation stop
                  if ((resolveMaterial && tLayer->hasMaterial()) || (resolveSensitive && tLayer->hasSensitive()) || tLayer == eLayer){
                      // taking the layer representation, not the approach surface, latter is resolved in the navigation stream
                      const Surface& tSurface = tLayer->surfaceRepresentation();
                      // intersect the layer @TODO should probably be surface on approach
                      Intersection lIntersection = tSurface.intersectionEstimate(gp,dir,true,bchk);
                      // (a) if the current layer is NOT the start layer - intersection is ok
                      if (tLayer != sLayer && lIntersection.valid){
                          lIntersections.push_back(LayerIntersection<T>(lIntersection,tLayer,&tSurface,0,pDir));
                          validPathLength = lIntersection.pathLength;
                      } else if (tLayer == sLayer) {
                          // (b) the current layer is the start layer - we need to cache it and check with the path length
                          //     this avoids potential punch-through to other side of
                          sLayerIntersection = lIntersection;
                          sLayerSurface      = &tSurface;
                      } else if (tLayer == eLayer) {
                          // (c) it is the end layer after all - provide it and break the loop
                          lIntersections.push_back(LayerIntersection<T>(lIntersection,tLayer,&tSurface,0,pDir));
                          break;
                      }
                  }
                  // move to next one or break because you reached the end layer
                  tLayer = (tLayer == eLayer ) ? nullptr :  tLayer->nextLayer(gp,dir);
              } while (tLayer);
          }

          // final check for compatibility of the start layer in order to avoid punch-through
          if (sLayer && sLayerIntersection.valid && sLayerIntersection.pathLength < validPathLength)
              lIntersections.push_back(LayerIntersection<T>(sLayerIntersection,sLayer,sLayerSurface,0,pDir));

      }
      // and the arbitraray layers
      if (m_confinedArbitraryLayers){
          // loop over the layers and intersect them
          for (auto& layer : (*m_confinedArbitraryLayers)){
              // intersections
              Intersection lIntersection = layer->surfaceRepresentation().intersectionEstimate(gp,dir,true,bchk);
              if (lIntersection.valid)
                  lIntersections.push_back(LayerIntersection<T>(lIntersection,layer.get(),&(layer->surfaceRepresentation()),0,pDir));
          }
      }

      // sort them accordingly to the path length
      std::sort(lIntersections.begin(),lIntersections.end());
      // and return
      return lIntersections;
  }

  /** Returns the boundary surfaces ordered in probability to hit them based on straight line intersection @TODO change hard-coded default */
  template <class T > std::vector< BoundaryIntersection<T> >
                      TrackingVolume::boundarySurfacesOrdered(const T& pars,
                                                              PropDirection pDir,
                                                              bool) const
  {

      // assign the direction
      const Vector3D dir = ( pDir == alongMomentum ? pars.momentum().unit() : Vector3D(-1*pars.momentum().unit()));
      // loop over boundarySurfaces and calculate the intersection
      std::vector< BoundaryIntersection<T> > bIntersections;
      auto& bSurfaces = boundarySurfaces();
      for (auto& bsIter : bSurfaces ){
          const BoundarySurface<TrackingVolume>* bSurface = bsIter.get();
          Intersection bsIntersection   = bSurface->surfaceRepresentation().intersectionEstimate(pars.position(),dir,true,false);
          if (bsIntersection.valid)
              bIntersections.push_back(BoundaryIntersection<T>(bsIntersection,bSurface,&(bSurface->surfaceRepresentation()),0,pDir));
      }
      // and now sort to get the closest
      std::sort(bIntersections.begin(),bIntersections.end());
      // and return
      return bIntersections;
  }

  inline GeometrySignature TrackingVolume::geometrySignature() const
  { return m_geometrySignature; }

  inline GeometryType TrackingVolume::geometryType() const
  { return m_geometryType; }

  inline void TrackingVolume::registerColorCode(unsigned int icolor) const
  { m_colorCode = icolor; }

  inline unsigned int TrackingVolume::colorCode() const
  { return m_colorCode; }

  inline const TrackingVolume* TrackingVolume::motherVolume() const
  { return m_motherVolume; }

  inline void TrackingVolume::setMotherVolume(const TrackingVolume* mvol) const
  { m_motherVolume = mvol; }

} // end of namespace

#endif // ACTS_DETECTOR_TRACKINGVOLUME_H

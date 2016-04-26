///////////////////////////////////////////////////////////////////
// TrackingGeometry.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_TRACKINGGEOMETRY_H
#define ACTS_DETECTOR_TRACKINGGEOMETRY_H

#include "ACTS/Utilities/GeometrySignature.h"
#include "ACTS/Utilities/Definitions.h"
#include <map>
#include <memory>

namespace Acts {
 
  class TrackingVolume;
  class DetachedTrackingVolume;
  class PerigeeSurface;
  class Layer;

  typedef std::shared_ptr< const TrackingVolume>         TrackingVolumePtr;
  typedef std::shared_ptr< const DetachedTrackingVolume> DetachedTrackingVolumePtr;
  typedef std::vector< DetachedTrackingVolumePtr >       DetachedVolumeVector;
   
  /** 
    @class TrackingGeometry
    
    The TrackingGeometry class is the owner of the constructed TrackingVolumes.
  
    It enables both, a global search for an asociatedVolume
    (respectively, if existing, a global search of an associated Layer or the next
    associated Layer), such as a continous navigation by BoundarySurfaces between 
    the confined TrackingVolumes.
    
    @author Andreas.Salzburger@cern.ch */
  
  class TrackingGeometry {

    /** Give the GeometryBuilder friend rights */  
    friend class TrackingGeometryBuilder;
  
    public :
      /** Constructor */
      TrackingGeometry(TrackingVolumePtr highestVolume);
      
      /** Destructor */
      ~TrackingGeometry();
      
      /** return the world */
      const TrackingVolume* highestTrackingVolume() const;
      
      /** return the lowest tracking Volume */
      const TrackingVolume* lowestTrackingVolume(const Vector3D& gp) const;

      /** return the vector of lowest detached tracking Volume(->overlaps) */
      const DetachedVolumeVector* lowestDetachedTrackingVolumes(const Vector3D& gp) const;

      /** return the lowest static tracking Volume */
      const TrackingVolume* lowestStaticTrackingVolume(const Vector3D& gp) const;
      
      /** return the tracking Volume by name, 0 if it doesn't exist */
      const TrackingVolume* trackingVolume(const std::string& name) const;
            
      /** Forward the associated Layer information */
      const Layer* associatedLayer(const Vector3D& gp) const;
      
      /** check position at volume boundary */
      bool atVolumeBoundary(const Vector3D& gp, const TrackingVolume* vol, double tol) const;
      
      /** check position at volume boundary + navigation link */
      bool atVolumeBoundary(const Vector3D& gp, const Vector3D& mom, const TrackingVolume* vol, 
			                const TrackingVolume*& nextVol, Acts::PropDirection dir, double tol) const;
      /** register the beam tube */
      void registerBeamTube(std::unique_ptr<const PerigeeSurface> beam) const;
      
    private:
      /** Geometry Builder busineess: the geometry builder has to sign*/
      void sign(GeometrySignature geosit, GeometryType geotype = Static) const;
      
      /** private method to register recursively the tracking volume & set the mother volume */
      void registerTrackingVolumes(const TrackingVolume& tvol, const TrackingVolume* mvol = nullptr, int lvl = 0);

      /** The known world - and the beamline */   
      TrackingVolumePtr  m_world;
      mutable std::unique_ptr<const PerigeeSurface>          m_beam;
      
      /** The Volumes in a map for later finding */
      std::map<const std::string, const TrackingVolume*>    m_trackingVolumes;

  };
  
} // end of namespace

#endif

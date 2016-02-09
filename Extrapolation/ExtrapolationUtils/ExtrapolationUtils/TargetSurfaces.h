///////////////////////////////////////////////////////////////////
// TransportJacobian.h, ATS project
///////////////////////////////////////////////////////////////////

#ifndef TRKEXUTILS_TARGETSURFACES_H
#define TRKEXUTILS_TARGETSURFACES_H


// CLHEP
#include "Algebra/AlgebraDefinitions.h"
#include "Surfaces/Surface.h"
#include "Surfaces/BoundaryCheck.h"
#include "ExtrapolationUtils/ExtrapolationCell.h"
#include "TrackParameters/TrackParameters.h"
#include "NeutralParameters/NeutralParameters.h"

class MsgStream;

namespace Ats {

  class TrackingVolume;

/** @class TargetSurfaces
  
   This class handles selection and distance estimate of surfaces whis can be reached
   in propagation/transport of particles.  

   @author sarka.todorova@cern.ch 
  
  */

  /** typedef for surface and volume navigation types */
  enum SurfNavigType
  {
    Target           = 0,
    BoundaryFrame    = 1,
    BoundaryDetached = 2,
    SensitiveLayer   = 3,     // sensitive layer
    MaterialLayer    = 4      // material layer
  };
  enum TVNavigType
  {
    Unknown          = 0,
    Frame            = 1,     // static volume  
    AlignableVolume  = 2,     // alignable volume
    DenseVolume      = 3      // dense alignable volume
  };

  /** target surface info ( navigation ) */
  struct TargetSurface{
    
    const Surface*         surf;        // surface
    BoundaryCheck          bcheck;      // boundary check
    SurfNavigType          sfType;      // surface type : boundary, layer ...
    unsigned int           index;       // layer or boundary index 
    const TrackingVolume*  assocVol;    // associated volume
    TVNavigType            volType;     // volume type : static, detached, dense ...
    double                 distanceAlongPath;     // latest calculated distance
    double                 distance;              // minimal distance estimate
    float                  signAbsDist;           // sign of minimal distance when available 
    Vector3D          intersection;          // global position 
    int                    status;          // navigation status (-1 crossed, 0 in loop, 1 ahead ) 
    
    /** Constructor */
    TargetSurface( const Surface* sf, BoundaryCheck bc, SurfNavigType stype, int ind,
		   const TrackingVolume* tVol, TVNavigType vtype ):
    surf(sf),bcheck(bc),sfType(stype),index(ind),assocVol(tVol),volType(vtype)
    {
      distanceAlongPath=0.; distance=0.; signAbsDist=0.; intersection=Vector3D(0.,0.,0.); status = 0; 
    }

    /** Destructor */
    ~TargetSurface(){};

    /** Distance info */
    void setDistance( double dAlongPath, double dMin, float sign) {
      distanceAlongPath = dAlongPath;
      distance          = dMin;
      signAbsDist       = sign;
    }

    /** Intersection info */
    void setPosition( Vector3D intPos ) {
      intersection = intPos;
    }

    /** Navigation status info */
    void setStatus( int st ) {
      status = st;
    }

    /** fast approximative update ( real distance > distance estimate ) */
    void fastDistanceUpdate( double step ) {
      // temporary check ( remove after development finished ) 
      // if ( fabs(distance) < step ) std::cout <<"fast distance update inappropriate " << std::endl;
      // end temporary

      distance = distance > 0 ? distance - step : distance + step;
      distanceAlongPath =  distanceAlongPath > 0 ?  distanceAlongPath - step : distanceAlongPath + step;

    }

  };

  typedef std::vector< TargetSurface > TargetSurfaceVector;
  typedef ExtrapolationCell<TrackParameters>   ExCellCharged;
  typedef ExtrapolationCell<NeutralParameters> ExCellNeutral;

  class TargetSurfaces  
    {
    public:
    
      /** Constructor */
      TargetSurfaces(){};
      
      /** Destructor */
      ~TargetSurfaces(){};

      /** Extract surfaces for charged propagation, step into new frame volume */
      Ats::ExtrapolationCode  setOnInput(Ats::ExCellCharged, const Ats::Surface* sf, BoundaryCheck bc) const;
      
      /** Extract surfaces for charged propagation, step into new frame volume */
      Ats::ExtrapolationCode  setOnInput(Vector3D position, Vector3D direction,
					 const Ats::TrackingVolume*, const Ats::Surface* sf, BoundaryCheck bc) const;

      /** Ordered intersections for neutral transport, step into new frame volume */
      TargetSurfaceVector  orderedIntersections(Ats::ExCellNeutral, const Ats::Surface* sf, BoundaryCheck bc) const;
      
      /** Ordered intersections for neutral transport, step into new frame volume */
      TargetSurfaceVector  orderedIntersections(Vector3D position, Vector3D direction,
						const Ats::TrackingVolume*, const Ats::Surface* sf, BoundaryCheck bc ) const;

      /** update of target surfaces at input or at frame volume boundary */
      bool initFrameVolume(Vector3D position, Vector3D direction,const Ats::TrackingVolume*) const;      

      /** intersections */
      void fillSolutions(int index, Vector3D gp, TargetSurfaceVector& solutions) const;

      /** distance reevaluation */
      bool checkDistance(Vector3D position, Vector3D direction, double nextStep) const;

      /** estimated distance along path to the nearest surface */
      double distanceToNext();

      /** step over intersection */
      bool flipDirection();
 
      /** index of nearest surface */
      int nextSf();

      /** number of intersections along path */
      unsigned numSf();

      /** debug mode */
      bool debugMode();

      /** set debug mode */
      void setDebugModeOn(){ m_debugMode = true;}

      /** set debug mode */
      void setDebugModeOff(){ m_debugMode = false;}

      /** current material volume/input for propagation */
      const Ats::TrackingVolume* currentDense();

      /** current material volume/input for propagation */
      const Ats::TrackingVolume* currentFrame();
      
    private:

      void evaluateInputDistance(Ats::TargetSurface& tt, Vector3D pos, Vector3D dir, bool base) const;
      void save(Ats::TargetSurface& tt, bool base) const;
      void findNext() const;
      TargetSurfaceVector  orderIntersections() const;
      bool updateDistance(int index, Ats::TargetSurface& tt, Vector3D position, Vector3D direction) const;

      mutable bool                              m_orderTrue;       // neutral(true)/charged(false)  
      mutable TargetSurfaceVector               m_baseSurfaces;    // surfaces to be followed all along the path;
      mutable std::vector<TargetSurfaceVector>  m_tempSurfaces;    // surfaces shadowed by an envelope  
      mutable TargetSurfaceVector               m_ordered;         // ordered intersections;

      mutable float                             m_tolerance;
      mutable Vector3D                     m_probePos;
      mutable Vector3D                     m_probeDir;
      mutable unsigned int                      m_numAlongPath;
      mutable int                               m_nextSf;
      mutable double                            m_distanceToNext;
      mutable double                            m_lastStep;
      mutable bool                              m_flipDirection;
      mutable double                            m_flipDistance;
  
      mutable const Ats::TrackingVolume*        m_currentFrame;
      mutable const Ats::TrackingVolume*        m_currentDense;

      mutable bool                              m_debugMode;        
      mutable bool                              m_absDist;
    };
           
  inline double TargetSurfaces::distanceToNext() { return m_distanceToNext; }

  inline bool TargetSurfaces::flipDirection() { return m_flipDirection; }

  inline int TargetSurfaces::nextSf() { return m_nextSf; }

  inline unsigned int TargetSurfaces::numSf() { return m_numAlongPath; }

  inline bool TargetSurfaces::debugMode() { return m_debugMode; }

  inline const Ats::TrackingVolume* TargetSurfaces::currentDense() { return m_currentDense; }
  
  inline const Ats::TrackingVolume* TargetSurfaces::currentFrame() { return m_currentFrame; }

} // end of namespace

#endif // TRKEXUTILS_TARGETSURFACES_H

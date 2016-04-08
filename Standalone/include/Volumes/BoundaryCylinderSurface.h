///////////////////////////////////////////////////////////////////
// BoundaryCylinderSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_BOUNDARYCYLIMDERSURFACE_H
#define ACTS_VOLUMES_BOUNDARYCYLIMDERSURFACE_H 1

// Geometry module
#include "Surfaces/CylinderSurface.h"
#include "Volumes/BoundarySurface.h"
// EventData module
#include "EventDataUtils/PropDirection.h"
// Gaudi
#include "GaudiKernel/SystemOfUnits.h"

namespace Acts {

    class Volume;

  /** 
   @class BoundaryCylinderSurface

   BoundaryCylinderSurface description inside the tracking realm,
   it extends the Surface description to make a surface being a boundary of a
   Acts::Volume
    
   @author Andreas.Salzburger@cern.ch
   */
      
  template <class T> class BoundaryCylinderSurface : 
                              virtual public BoundarySurface<T>, public CylinderSurface {
  
    typedef std::shared_ptr<const T> VolumePtr; 
    typedef BinnedArray< VolumePtr > VolumeArray; 
  
    public:
     /** Default Constructor - needed for pool and inherited classes */
     BoundaryCylinderSurface():
      BoundarySurface<T>(),
      CylinderSurface()
     {}
     
     /** Copy constructor */                            
     BoundaryCylinderSurface(const BoundaryCylinderSurface<T>& bcs) :
       BoundarySurface<T>(bcs),
       CylinderSurface(bcs)
     {}
     
     /** Constructor for a Boundary with exact two Volumes attached to it*/
     BoundaryCylinderSurface(const T* inside, const T* outside, const CylinderSurface& csf) :
       BoundarySurface<T>(inside, outside),
       CylinderSurface(csf)
     {}
       
     /** Constructor for a Boundary with two VolumeArrays attached to it*/
     BoundaryCylinderSurface(std::shared_ptr<const VolumeArray> insideArray, std::shared_ptr<const VolumeArray> outsideArray, const CylinderSurface& csf) :
       BoundarySurface<T>(insideArray, outsideArray),
       CylinderSurface(csf)
     {}       
          
     /** Copy constructor with a shift */
     BoundaryCylinderSurface(const T* inside, const T* outside, const CylinderSurface& csf, const Transform3D& tr) :
       BoundarySurface<T>(inside,outside),
       CylinderSurface(csf,tr)
     {}
          
     /**Virtual Destructor*/
     virtual ~BoundaryCylinderSurface(){}
     
     /** Get the next Volume depending on position, momentum, dir
      of the TrackParameters and the requested direction */
     const T* attachedVolume(const Vector3D& pos, const Vector3D& mom, PropDirection dir) const override;
                                          
     /** The Surface Representation of this */
     const Surface& surfaceRepresentation() const override;
     
     /**Assignment operator - assignment of boundary surfaces is forbidden */
     BoundaryCylinderSurface& operator=(const BoundaryCylinderSurface& vol) = delete;
     
   protected:
                             
  };

  template <class T> inline const Surface& BoundaryCylinderSurface<T>::surfaceRepresentation() const { return *this; }
  
  template <class T> inline const T* BoundaryCylinderSurface<T>::attachedVolume(const Vector3D& pos,
                                                                                const Vector3D& mom,
                                                                                PropDirection dir) const
  {
    const T* attVolume = nullptr;
    // it's a cylinder surface, we need to get the local normal vector
    // -> bring the position vector into the cylinder frame (if needed)
    Vector3D normalSf = Surface::m_transform ? Surface::m_transform->inverse()*pos : pos;
    normalSf[2] = 0.;
    // ->bring the momentum vector into the cylinder frame (if needed)
    Vector3D momentumSf = Surface::m_transform ? Surface::m_transform->inverse()*mom : mom;
    // dot product to decide between inside and outside)
    if ( normalSf.unit().dot(dir*momentumSf) > 0. ){
        attVolume = BoundarySurface<T>::m_outsideVolume;
        if (BoundarySurface<T>::m_outsideVolumeArray.get())
          attVolume = BoundarySurface<T>::m_outsideVolumeArray->object(pos).get();    
    } else {
        attVolume = BoundarySurface<T>::m_insideVolume;
        if (BoundarySurface<T>::m_insideVolumeArray.get())
          attVolume = BoundarySurface<T>::m_insideVolumeArray->object(pos).get();
    }
    // return what you have
    return attVolume;  
  }
  
} // end of namespace Acts

#endif // ACTS_VOLUMES_BOUNDARYCYLIMDERSURFACE_H


///////////////////////////////////////////////////////////////////
// BoundaryDiscSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_BOUNDARYDISCSURFACE_H
#define ACTS_VOLUMES_BOUNDARYDISCSURFACE_H 1

// Geometry module
#include "ACTS/Utilities/Definitions.h"
#include "ACTS/Surfaces/DiscSurface.h"
#include "ACTS/Volumes/BoundarySurface.h"
// EventData module
#include "ACTS/Utilities/Definitions.h"
// Core module

namespace Acts {

   class Volume;

  /**
    @class BoundaryDiscSurface

    BoundaryDiscSurface description inside the tracking realm,
    it extends the DiscSurface description to make a surface being a boundary of a
    Acts::Volume (used for cylindrical shape).
    It inherits from BoundarySurface to get the interface of boundaries.

    @author Andreas.Salzburger@cern.ch
   */

  template <class T> class BoundaryDiscSurface :
                       virtual public BoundarySurface<T>, public DiscSurface {

    typedef std::shared_ptr<const T> VolumePtr;
    typedef BinnedArray< VolumePtr > VolumeArray;

    public:
     /** Default Constructor - needed for pool and inherited classes */
     BoundaryDiscSurface():
      BoundarySurface<T>(),
      DiscSurface()
     {}

     /** Copy constructor */
     BoundaryDiscSurface(const BoundaryDiscSurface<T>& bds) :
       BoundarySurface<T>(bds),
       DiscSurface(bds)
     {}

     /** Constructor for a Boundary with exact two Volumes attached to it*/
     BoundaryDiscSurface(const T* inside, const T* outside, const DiscSurface& dsf) :
       BoundarySurface<T>(inside, outside),
       DiscSurface(dsf)
     {}

     /** Constructor for a Boundary with two VolumeArrays attached to it*/
     BoundaryDiscSurface(std::shared_ptr<const VolumeArray> insideArray, std::shared_ptr<const VolumeArray> outsideArray, const DiscSurface& dsf) :
       BoundarySurface<T>(insideArray, outsideArray),
       DiscSurface(dsf)
     {}

     /** Copy constructor with a shift */
     BoundaryDiscSurface(const T* inside, const T* outside, const DiscSurface& dsf, const Transform3D& tr) :
       BoundarySurface<T>(inside,outside),
       DiscSurface(dsf,tr)
     {}

     /** The Surface Representation of this */
     const Surface& surfaceRepresentation() const override;

     /**Virtual Destructor*/
     virtual ~BoundaryDiscSurface(){}

     /**Assignment operator - assignment of boundary surfaces is forbidden */
     BoundaryDiscSurface& operator=(const BoundaryDiscSurface& vol) = delete;

   protected:

  };

template <class T> inline const Surface& BoundaryDiscSurface<T>::surfaceRepresentation() const { return *this; }

} // end of namespace Acts

#endif // ACTS_VOLUMES_BOUNDARYDISCSURFACE_H


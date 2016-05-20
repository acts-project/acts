// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BoundarySubtractedCylinderSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_BOUNDARYSUBTRACTEDCYLIMDERSURFACE_H
#define ACTS_VOLUMES_BOUNDARYSUBTRACTEDCYLIMDERSURFACE_H 1

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Surfaces/SubtractedCylinderSurface.hpp"
#include "ACTS/Volumes/BoundarySurface.hpp"

namespace Acts {

  class Volume;

  /**
   @class BoundarySubtractedCylinderSurface

   BoundarySubtractedCylinderSurface description inside the tracking realm,
   it extends the Surface description to make a surface being a boundary of a
   Acts::Volume

   */

  template <class T> class BoundarySubtractedCylinderSurface :
                              virtual public BoundarySurface<T>, public SubtractedCylinderSurface {

    /** typedef the BinnedArray */
    typedef BinnedArray<T> VolumeArray;

    public:
     /** Default Constructor - needed for pool and inherited classes */
     BoundarySubtractedCylinderSurface():
      BoundarySurface<T>(),
      SubtractedCylinderSurface()
     {}

     /** Copy constructor */
     BoundarySubtractedCylinderSurface(const BoundarySubtractedCylinderSurface<T>& bcs) :
       BoundarySurface<T>(bcs),
       SubtractedCylinderSurface(bcs)
     {}

     /** Constructor for a Boundary with exact two Volumes attached to it*/
     BoundarySubtractedCylinderSurface(const T* inside, const T* outside, const SubtractedCylinderSurface& csf) :
       BoundarySurface<T>(inside, outside),
       SubtractedCylinderSurface(csf)
     {}

     /** Constructor for a Boundary with two VolumeArrays attached to it*/
     BoundarySubtractedCylinderSurface(std::shared_ptr<VolumeArray> insideArray, std::shared_ptr<VolumeArray> outsideArray, const SubtractedCylinderSurface& csf) :
       BoundarySurface<T>(insideArray, outsideArray),
       SubtractedCylinderSurface(csf)
     {}

     /** Copy constructor with a shift */
     BoundarySubtractedCylinderSurface(const T* inside, const T* outside, const SubtractedCylinderSurface& csf, const Transform3D& tr) :
       BoundarySurface<T>(inside,outside),
       SubtractedCylinderSurface(csf,tr)
     {}

     /**Virtual Destructor*/
     virtual ~BoundarySubtractedCylinderSurface()
     {}

     /** Get the next Volume depending on GlobalPosition, GlobalMomentum, dir on the TrackParameters and the requested direction */
     const T* attachedVolume(const Vector3D& pos, const Vector3D& mom, PropDirection dir) const override;

     /** The Surface Representation of this */
     const Surface& surfaceRepresentation() const override;

     /**Assignment operator - forbidden */
     BoundarySubtractedCylinderSurface& operator=(const BoundarySubtractedCylinderSurface& vol) = delete ;

  };

  template <class T> inline const Surface& BoundarySubtractedCylinderSurface<T>::surfaceRepresentation() const { return *this; }

  template <class T> inline const T* BoundarySubtractedCylinderSurface<T>::attachedVolume(const Vector3D& pos,
                                                                                          const Vector3D& mom,
                                                                                          PropDirection dir) const
  {
      const T* attVolume = 0;
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

#endif // ACTS_VOLUMES_BOUNDARYSUBTRACTEDCYLINDERSURFACE_H

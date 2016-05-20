// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BoundarySurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_BOUNDARYSURFACE_H
#define ACTS_VOLUMES_BOUNDARYSURFACE_H

#include "ACTS/Volumes/Volume.hpp"
#include "ACTS/Volumes/BoundarySurfaceFace.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/BinnedArray0D.hpp"
#include "ACTS/Utilities/Definitions.hpp"


namespace Acts {

  class Surface;

  /**
   @class BoundarySurface

   Description of a BoundarySurface for volumes in the Tracking geometry.
   It extends the Surface description to make a surface being a boundary of a
   Acts::AbstractVolume & Acts::TrackingVolume.

   To avoid dynamic_cast operations the BoundarySurface class is realized as a templated class,
   with the Volume type as the template argument.

   A Acts::BoundarySurface can have an inside Volume and an outside Volume, resp.
   a Acts::BinnedArray for inside or outside direction.

   Inside and outside are by this referring to the normal vector definition of the Surface,
   and not necessarily corresponding to the Volume inside/outside definition.

   The GeometryBuilder as defined in the GeometryTools Package is declared
   to be friend, so that it can glue Volumes together by sharing the same
   Boundary Surface.

   @author Andreas.Salzburger@cern.ch
  */

  template <class T> class BoundarySurface {

    /** delcare the TrackingVolume as friend */
    friend T;

    typedef std::shared_ptr<const T> VolumePtr;
    typedef BinnedArray< VolumePtr > VolumeArray;

    public:
     /** Default Constructor - needed for pool and inherited classes */
     BoundarySurface() :
       m_insideVolume(nullptr),
       m_outsideVolume(nullptr),
       m_insideVolumeArray(nullptr),
       m_outsideVolumeArray(nullptr)
     {}

     /** Constructor for a Boundary with exact two Volumes attached to it - usually used in a volume constructor*/
     BoundarySurface(const T* inside, const T* outside) :
       m_insideVolume(inside),
       m_outsideVolume(outside),
       m_insideVolumeArray(nullptr),
       m_outsideVolumeArray(nullptr)
     {}

     /** Constructor for a Boundary with exact two Volumes attached to it - for already constructed volumes */
     BoundarySurface(VolumePtr inside, VolumePtr outside) :
       m_insideVolume(nullptr),
       m_outsideVolume(nullptr),
       m_insideVolumeArray(new BinnedArray0D<VolumePtr>(inside)),
       m_outsideVolumeArray(new BinnedArray0D<VolumePtr>(outside))
     {}

     /** Constructor for a Boundary with exact two Volumes attached to it*/
     BoundarySurface(std::shared_ptr<const VolumeArray> insideArray, std::shared_ptr<const VolumeArray> outsideArray) :
       m_insideVolume(nullptr),
       m_outsideVolume(nullptr),
       m_insideVolumeArray(insideArray),
       m_outsideVolumeArray(outsideArray)
     {}

     /** Get the next Volume depending on GlobalPosition, GlobalMomentum, dir on the TrackParameters and the requested direction
         - the position, momentum are assumed to be on surface  */
     virtual const T* attachedVolume(const Vector3D& pos, const Vector3D& mom, PropDirection dir) const;

     /** templated onBoundary method */
     template <class P> bool onBoundary(const P& pars) const
     { return surfaceRepresentation().onSurface(pars); }

     /** The Surface Representation of this */
     virtual const Surface& surfaceRepresentation() const = 0;

     /**Virtual Destructor*/
     virtual ~BoundarySurface(){}

   protected:
     /** attach a Volume to this BoundarySurface - will always be done as an 0D array
         - the owning volume can call this as it is fiend */
     void attachVolume(VolumePtr volume, BoundaryOrientation inout) const;

     /** attach a  VolumeArray to this BoundarySurface
         - the owing volume can call this as it is friend */
     void attachVolumeArray(std::shared_ptr<const VolumeArray> volumes, BoundaryOrientation inout) const;

     mutable const T*                            m_insideVolume;
     mutable const T*                            m_outsideVolume;
     mutable std::shared_ptr<const VolumeArray>  m_insideVolumeArray;
     mutable std::shared_ptr<const VolumeArray>  m_outsideVolumeArray;

  };

  template <class T> void BoundarySurface<T>::attachVolume(VolumePtr volume, BoundaryOrientation inout) const {
      if (inout == insideVolume) m_insideVolumeArray = std::shared_ptr<const VolumeArray>(new BinnedArray0D<VolumePtr>(volume));
      else m_outsideVolumeArray = std::shared_ptr<const VolumeArray>(new BinnedArray0D<VolumePtr>(volume));
  }

  template <class T> void BoundarySurface<T>::attachVolumeArray(const std::shared_ptr<const VolumeArray> volumes, BoundaryOrientation inout) const {
      if (inout == insideVolume) m_insideVolumeArray = volumes;
      else m_outsideVolumeArray = volumes;
  }

  template <class T> const T* BoundarySurface<T>::attachedVolume(const Vector3D& pos, const Vector3D& mom, PropDirection dir) const
  {

    const T* attVolume = nullptr;
    // dot product with normal vector to distinguish inside/outside
    if ( (surfaceRepresentation().normal(pos)).dot(dir*mom) > 0.)
        attVolume = m_outsideVolumeArray ? m_outsideVolumeArray->object(pos).get() : m_outsideVolume;
    else
        attVolume = m_insideVolumeArray ? m_insideVolumeArray->object(pos).get() : m_insideVolume;
   return attVolume;
  }
} // end of namespace Acts

#endif // ACTS_VOLUMES_BOUNDARYSURFACE_H

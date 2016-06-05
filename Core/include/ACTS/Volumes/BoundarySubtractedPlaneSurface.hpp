// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BoundarySubtractedPlaneSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_BOUNDARYSUBTRACTEDPLANESURFACE_H
#define ACTS_VOLUMES_BOUNDARYSUBTRACTEDPLANESURFACE_H 1

#include "ACTS/Surfaces/SubtractedPlaneSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Volumes/BoundarySurface.hpp"

namespace Acts {

class Volume;

/**
 @class BoundarySubtractedPlaneSurface

 BoundarySubtractedPlaneSurface description inside the tracking realm,
 it extends the SubtractedPlaneSurface description to make a surface being a
 boundary of a
 Acts::Volume (used for all volume shapes).
 It inherits from BoundarySurface to get the interface of boundaries.

*/

template <class T>
class BoundarySubtractedPlaneSurface : virtual public BoundarySurface<T>,
                                       public SubtractedPlaneSurface
{
  /** typedef the BinnedArray */
  typedef BinnedArray<T> VolumeArray;

public:
  /** Default Constructor - needed for pool and inherited classes */
  BoundarySubtractedPlaneSurface()
    : BoundarySurface<T>(), SubtractedPlaneSurface()
  {
  }

  /** Copy constructor */
  BoundarySubtractedPlaneSurface(const BoundarySubtractedPlaneSurface<T>& bps)
    : BoundarySurface<T>(bps), SubtractedPlaneSurface(bps)
  {
  }

  /** Constructor for a Boundary with exact two Volumes attached to it*/
  BoundarySubtractedPlaneSurface(const T*                      inside,
                                 const T*                      outside,
                                 const SubtractedPlaneSurface& psf)
    : BoundarySurface<T>(inside, outside), SubtractedPlaneSurface(psf)
  {
  }

  /** Constructor for a Boundary with two VolumeArrays attached to it*/
  BoundarySubtractedPlaneSurface(std::shared_ptr<VolumeArray>  insideArray,
                                 std::shared_ptr<VolumeArray>  outsideArray,
                                 const SubtractedPlaneSurface& psf)
    : BoundarySurface<T>(insideArray, outsideArray), SubtractedPlaneSurface(psf)
  {
  }

  /** Copy constructor with a shift */
  BoundarySubtractedPlaneSurface(const T*                      inside,
                                 const T*                      outside,
                                 const SubtractedPlaneSurface& psf,
                                 const Transform3D&            tr)
    : BoundarySurface<T>(inside, outside), SubtractedPlaneSurface(psf, tr)
  {
  }

  /** The Surface Representation of this */
  const Surface&
  surfaceRepresentation() const override;

  /**Virtual Destructor*/
  virtual ~BoundarySubtractedPlaneSurface() {}
  /**Assignment operator - assignment is forbidden */
  BoundarySubtractedPlaneSurface&
  operator=(const BoundarySubtractedPlaneSurface& vol)
      = delete;
};

template <class T>
inline const Surface&
BoundarySubtractedPlaneSurface<T>::surfaceRepresentation() const
{
  return *this;
}

}  // end of namespace Acts

#endif  // ACTS_VOLUMES_BOUNDARYSUBTRACTEDPLANESURFACE_H

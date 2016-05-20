// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AbstractVolume.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_ABSTRACTVOLUME_H
#define ACTS_VOLUMES_ABSTRACTVOLUME_H

// Geometry module
#include "ACTS/Volumes/Volume.hpp"
#include "ACTS/Volumes/BoundarySurface.hpp"
// Core module
#include <vector>
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

  class VolumeBounds;
  typedef std::shared_ptr<const VolumeBounds> VolumeBoundsPtr;

  /**
   @class AbstractVolume

   AbstractVolume description inside the tracking realm. This is the purely geometrical object volume.

   The Acts::AbstractVolume is constructed by giving a pointer to a Transform3D
   and a pointer to Acts::VolumeBounds, this implies that the ownership of the
   objects pointed to is passed as well. For memory optimisation, the AbstractVolume can also be
   constructed with shared_ptr objects.

   A Acts::AbstractVolume is at first a collection class of Acts::BoundarySurface,
   the vector of Acts::BoundarySurface is returned by the Acts::VolumeBounds that
   carry a decompose method.

   Boundary surfaces can be shared between AbstractVolumes to enhance automatic navigation
   between AbstractVolumes, therefor they are reference counted by a std::shared_ptr holder class.

   @image html VolumeShapes.gif

   @author Andreas.Salzburger@cern.ch
   */

  class AbstractVolume : public Volume {

    public:
     /**Constructor with Transform3D*, VolumeBounds*, passing ownership */
     AbstractVolume(Transform3D* htrans, const VolumeBounds* volbounds);

     /**Constructor with shared Transform3D*, VolumeBounds* */
     AbstractVolume(std::shared_ptr<Transform3D> htrans, VolumeBoundsPtr volbounds);

     /**Copy constructor - deleted */
     AbstractVolume(const AbstractVolume& vol) = delete;

     /**Virtual Destructor*/
     virtual ~AbstractVolume();

     /**Assignment operator*/
     AbstractVolume& operator=(const AbstractVolume& vol);

     /** Method to return the BoundarySurfaces */
     const std::vector< std::shared_ptr<const BoundarySurface<AbstractVolume> > >& boundarySurfaces() const;

   private:
     /**Default Constructor - needed for pool and inherited classes */
     AbstractVolume();

     /**Private method to create BoundarySurfaces */
     void createBoundarySurfaces();

     std::vector< std::shared_ptr<const BoundarySurface<AbstractVolume> > >* m_boundarySurfaces;  //!< boundary Surfaces

  };

} // end of namespace

#endif // ACTS_VOLUMES_ABSTRACTVOLUME_H


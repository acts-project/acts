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

#include "ACTS/Volumes/BoundarySurface.hpp"
#include "ACTS/Volumes/Volume.hpp"
#include <vector>
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

class VolumeBounds;
typedef std::shared_ptr<const VolumeBounds> VolumeBoundsPtr;

typedef std::shared_ptr<const BoundarySurface<AbstractVolume> > BoundarySurfacePtr;
/// @class AbstractVolume
///
/// AbstractVolume description inside the tracking realm. This is the purely
/// geometrical object volume.
///
/// The Acts::AbstractVolume is constructed by giving a pointer to a Transform3D
/// and a pointer to Acts::VolumeBounds, this implies that the ownership of the
/// objects pointed to is passed as well. For memory optimisation, the
/// AbstractVolume can also be
/// constructed with shared_ptr objects.
///
/// A Acts::AbstractVolume is at first a collection class of Acts::BoundarySurface,
/// the vector of Acts::BoundarySurface is returned by the Acts::VolumeBounds that
/// carry a decompose method.
///
/// Boundary surfaces can be shared between AbstractVolumes to enhance automatic
/// navigation
/// between AbstractVolumes, therefor they are reference counted by a
/// std::shared_ptr holder class.
///
/// @image html VolumeShapes.gif


class AbstractVolume : public Volume
{
    
public:
  /// Constructor with shared Transform3D*, VolumeBounds* 
  /// @param htrans is the transform 3D the positions the volume in global frame
  /// @param volbounds is the boundary definition
  AbstractVolume(std::shared_ptr<Transform3D> htrans,
                 VolumeBoundsPtr              volbounds);

  /// Copy constructor - deleted 
  AbstractVolume(const AbstractVolume& vol) = delete;

  /// Default Constructor - deleted
  AbstractVolume() = delete;

  // Virtual Destructor
  virtual ~AbstractVolume();

  /// Assignment operator - deleted
  AbstractVolume&
  operator=(const AbstractVolume& vol) = delete;

  /// Method to return the BoundarySurfaces 
  const std::vector< BoundarySurfacePtr >& boundarySurfaces() const;

private:

  /// Private method to create BoundarySurfaces 
  void
  createBoundarySurfaces();

  /// boundary Surfaces for this volume
  std::vector< BoundarySurfacePtr > m_boundarySurfaces;  
};

}  // end of namespace

#endif  // ACTS_VOLUMES_ABSTRACTVOLUME_H

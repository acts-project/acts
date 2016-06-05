// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GlueVolumesDescriptor.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_GLUEVOLUMESDESCRIPTOR_H
#define ACTS_DETECTOR_GLUEVOLUMESDESCRIPTOR_H 1

// Geometry module
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Volumes/BoundarySurfaceFace.hpp"
// STL
#include <map>
#include <memory>
#include <vector>

namespace Acts {

class TrackingVolume;
typedef std::shared_ptr<const TrackingVolume> TrackingVolumePtr;
typedef BinnedArray<TrackingVolumePtr>        TrackingVolumeArray;

/** @class GlueVolumesDescriptor

   Descriptor class to hold GlueVolumes of a TrackingGeometry object.
   Should ease the wrapping of a TrackingGeometry object describing one Detector
   by another one.

   */

class GlueVolumesDescriptor
{
public:
  /**Constructor */
  GlueVolumesDescriptor() {}
  /**Constructor - with arguments*/
  GlueVolumesDescriptor(
      const std::map<BoundarySurfaceFace,
                     std::shared_ptr<const TrackingVolumeArray>>& gvs);

  /**Desctructor */
  ~GlueVolumesDescriptor() {}
  /** register the volumes */
  void
  registerGlueVolumes(Acts::BoundarySurfaceFace                  bsf,
                      std::shared_ptr<const TrackingVolumeArray> gvs) const;

  /** retrieve them again */
  std::shared_ptr<const TrackingVolumeArray>
      glueVolumes(BoundarySurfaceFace) const;

  /** retrieve the available Glue Faces */
  const std::vector<BoundarySurfaceFace>&
  glueFaces() const;

  /** dump it to the screen */
  std::string
  screenOutput() const;

private:
  mutable std::map<BoundarySurfaceFace,
                   std::shared_ptr<const TrackingVolumeArray>>
                                           m_glueVolumes;
  mutable std::vector<BoundarySurfaceFace> m_glueFaces;
};

inline const std::vector<BoundarySurfaceFace>&
GlueVolumesDescriptor::glueFaces() const
{
  return m_glueFaces;
}

/**Overload of << operator for std::ostream for debug output*/
std::ostream&
operator<<(std::ostream& sl, const GlueVolumesDescriptor& mprop);
}

#endif

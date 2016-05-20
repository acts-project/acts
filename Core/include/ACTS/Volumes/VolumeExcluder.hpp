// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// VolumeExcluder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_VOLUMEEXCLUDER_H
#define ACTS_VOLUMES_VOLUMEEXCLUDER_H 1

// Geometry module
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/AreaExcluder.hpp"
#include "ACTS/Volumes/Volume.hpp"
// Core module

namespace Acts {

  class AreaExcluder;

/** @class VolumeExcluder
    removes explicit dependence of Subtracted*Surface on Volumes

   @author sarka.todorova@cern.ch
  */

   class VolumeExcluder: public AreaExcluder {

      public:
        /** Default constructor */
        VolumeExcluder();

        /** Explicit constructor with volume */
        VolumeExcluder(Volume* vol);

        /** copy constructor */
        VolumeExcluder(const VolumeExcluder& ex);

        /** Destructor */
        virtual ~VolumeExcluder();

        /** Assignment operator */
        VolumeExcluder& operator=(const VolumeExcluder &vol);

        /** Pseudo-constructor */
        VolumeExcluder* clone() const;

        /** First bin from global position */
        bool inside(const Vector3D& gp, double tol=0.) const;

        /** Acces the subtracted volume */
        Volume* volume() const;

        /** Output Method for std::ostream, to be overloaded by child classes */
        std::ostream& dump(std::ostream& sl) const;

     private:
        Volume* m_vol;

   };

   inline bool VolumeExcluder::inside(const Vector3D& gp, double tol) const
    {  return ( m_vol->inside(gp,tol) ); }

   inline Volume* VolumeExcluder::volume() const
    {  return ( m_vol ); }

} // end of namespace Acts

#endif // ACTS_VOLUMES_VOLUMEEXCLUDER


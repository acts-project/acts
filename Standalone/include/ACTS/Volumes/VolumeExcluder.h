///////////////////////////////////////////////////////////////////
// VolumeExcluder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_VOLUMEEXCLUDER_H
#define ACTS_VOLUMES_VOLUMEEXCLUDER_H 1

// Geometry module
#include "ACTS/GeometryUtils/AreaExcluder.h"
#include "ACTS/Volumes/Volume.h"
// Core module
#include "ACTS/Utilities/AlgebraDefinitions.h"

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


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
#ifndef ACTS_VOLUMES_VOLUME_H
#define ACTS_VOLUMES_VOLUME_H 1

// STL include(s)
#include <memory>

#include "ACTS/Utilities/Definitions.hpp"
// Geometry module
#include "ACTS/Utilities/GeometryStatics.hpp"
#include "ACTS/Utilities/GeometryObject.hpp"
// Core module

namespace Acts {

  class VolumeBounds;
  typedef std::shared_ptr<const VolumeBounds> VolumeBoundsPtr;

  /** @class Volume

    It inhertis of GeometryObject for TDD identification

    Base class for all volumes inside the tracking realm, it defines
    the interface for inherited Volume classes
    regarding the geometrical information.

    @author Andreas.Salzburger@cern.ch
    */

    class Volume : public virtual GeometryObject {

      public:
        /** Default constructor */
        Volume();

        /** Expizit constructor with arguments */
        Volume(Transform3D* htrans, const VolumeBounds* volBounds);

        /** Expizit constructor with shared arguments */
        Volume(std::shared_ptr<Transform3D> htrans, VolumeBoundsPtr volBounds);

        /** Copy Constructor */
        Volume(const Volume& vol);

        /** Copy Constructor */
        Volume(const Volume& vol, const Transform3D& shift);

        /** Destructor */
        virtual ~Volume();

        /** Assignment operator */
        Volume& operator=(const Volume& vol);

        /** Pseudo-constructor */
        virtual Volume* clone() const;

        /** Return methods for geometry transform */
        const Transform3D& transform() const;

        /** returns the center of the volume */
        const Vector3D& center() const;

        /** returns the volumeBounds() */
        const VolumeBounds& volumeBounds() const;

        /** Inside() method for checks */
        bool inside(const Vector3D& gp, double tol=0.) const;

        /** The binning position method - as default the center is given, but may be overloaded */
        virtual Vector3D binningPosition(BinningValue bValue) const override;

      protected:
        std::shared_ptr<Transform3D>             m_transform;         //!< Transform3D
        mutable Vector3D*                        m_center;            //!< center position of the surface
        VolumeBoundsPtr                          m_volumeBounds;      //!< the volumeBounds
    };

    inline const Transform3D& Volume::transform() const
    {  if (m_transform.get()) return(*(m_transform.get()));
       return Acts::s_idTransform;
    }

    inline const Vector3D& Volume::center() const
    {
     if (m_center) return (*m_center);
     if (!m_center && m_transform.get()){
        m_center = new Vector3D(m_transform->translation());
        return(*m_center);
      }
      return Acts::s_origin;
    }

    inline const VolumeBounds& Volume::volumeBounds() const
    {  return (*(m_volumeBounds.get())); }


/**Overload of << operator for std::ostream for debug output*/
std::ostream& operator << ( std::ostream& sl, const Volume& vol);


} // end of namespace Acts

#endif // ACTS_VOLUMES_VOLUME_H

// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AreaExcluder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_AREAEXCLUDER_H
#define ACTS_GEOMETRYUTILS_AREAEXCLUDER_H 1

// Core module
#include "Definitions.hpp"

namespace Acts {

/** @class AreaExcluder
    Pure abstract base class

  */

   class AreaExcluder {

      public:
        /** Default constructor */
        AreaExcluder(){}

        /** Virtual Destructor */
        virtual ~AreaExcluder(){}

        /** Implizit Constructor */
        virtual AreaExcluder* clone() const = 0;

        /** First bin from global position */
        virtual bool inside(const Vector3D& gp, double tol=0.) const = 0;

   };

} // end of namespace Acts

#endif // TRKDETDESCRUTILS_AREAEXCLUDER_H


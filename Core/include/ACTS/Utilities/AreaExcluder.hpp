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

    @author sarka.todorova@cern.ch
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


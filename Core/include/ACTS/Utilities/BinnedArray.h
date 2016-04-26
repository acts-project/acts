///////////////////////////////////////////////////////////////////
// BinnedArray.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_BINNEDARRAY_H
#define ACTS_GEOMETRYUTILS_BINNEDARRAY_H 1

// Geometry module
#include "ACTS/Utilities/BinUtility.h"
// Core moudle
#include <vector>
#include "Definitions.h"

namespace Acts {

/** @class BinnedArray

   Pure virtual base class for Binned Array to avoid map searches
   - there is only one restriction:
     T must be of pointer type in order to be initialized withh nullptr
     and to allow for nulpptr return type

   @author Andreas.Salzburger@cern.ch
   */

  template <class T> class BinnedArray {

    public:
      /** Default Constructor - needed for inherited classes */
      BinnedArray() {}

      /** Default Constructor - needed for inherited classes */
      BinnedArray(const BinnedArray<T>& ba) :
        m_arrayObjects(ba.m_arrayObjects)
      {}

      /**Virtual Destructor*/
      virtual ~BinnedArray(){}

      /** Implicit constructor */
      virtual BinnedArray* clone() const = 0;

      /** Returns the object in the associated bin according the local position  */
      virtual T object(const Vector2D& lp) const = 0;

      /** Returns the object in the associated bin according the global position  */
      virtual T object(const Vector3D& gp) const = 0;

      /** Returns the pointer to the templated class object from the BinnedArray - entry point*/
      virtual T entryObject(const Vector3D&) const = 0;
      
      /** Return all objects of the Array - in one */
      virtual const std::vector< T >& arrayObjects() const { return m_arrayObjects; };

      /** Return the BinUtility*/
      virtual const BinUtility* binUtility() const = 0;

    protected:
      std::vector< T >    m_arrayObjects;

  };

} // end of namespace Acts

#endif // ACTS_GEOMETRYUTILS_BINNEDARRAY_H

// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BinnedArray0D.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_BINNEDARRAY0D_H
#define ACTS_GEOMETRYUTILS_BINNEDARRAY0D_H 1

// Geometry module
#include "ACTS/Utilities/BinUtility.hpp"
// STL
#include <vector>

class MsgStream;

namespace Acts {

/** @class BinnedArray0D

   BinnedArray0D - BinnedArray flavour with just one entry,
   allows to have the same structure for just one single
   contained object.

   @author Andreas.Salzburger@cern.ch
   */

  template <class T> class BinnedArray0D : public BinnedArray<T> {

    public:
     /**Default Constructor - needed for inherited classes */
     BinnedArray0D():
       BinnedArray<T>()
     {}

     /**Constructor  from oneobject  */
     BinnedArray0D(T obj):
       BinnedArray<T>(),
       m_object(obj)
     {
       BinnedArray<T>::m_arrayObjects.push_back(obj);
     }

     /**Copy Constructor - copies only pointers !*/
     BinnedArray0D(const BinnedArray0D& barr):
       BinnedArray<T>(barr),
       m_object(barr.m_object)
      {}

     /**Assignment operator*/
     BinnedArray0D& operator=(const BinnedArray0D& barr)
     {
       if (this != &barr){
          m_object       = barr.m_object;
       }
        return *this;
     }

     /** Implizit Constructor */
     BinnedArray0D* clone() const final
     { return new BinnedArray0D(*this); }

     /**Virtual Destructor*/
     ~BinnedArray0D(){}

     /** Returns the single contained object */
     T object(const Vector2D&) const final { return m_object; }

     /** Returns the single contained object */
     T object(const Vector3D&) const final { return m_object; }

     /** Returns the single contained object */
     T entryObject(const Vector3D&) const final { return m_object; }

     /** Return the BinUtility*/
     const BinUtility* binUtility() const { return nullptr; }

    private:
     T          m_object;       //!< the single contained object

  };

} // end of namespace Acts

#endif // ACTS_GEOMETRYUTILS_BINNEDARRAY0D_H

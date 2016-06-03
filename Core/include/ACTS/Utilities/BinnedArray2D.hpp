// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BinnedArray2D.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_BINNEDARRAY2D_H
#define ACTS_GEOMETRYUTILS_BINNEDARRAY2D_H 1

// Geometry module
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
// STL
#include <vector>

class MsgStream;

namespace Acts {
    
    /** @class BinnedArray2D
     
     Avoiding a map search, the templated BinnedArray class can help
     ordereing geometrical objects by providing a dedicated BinUtility.
     
     dedicated for 2-dim regular binning
     
     */
    
    template <class T> class BinnedArray2D : public BinnedArray<T> {
        
    public:
        
        /**Default Constructor - needed for inherited classes */
        BinnedArray2D() :
         BinnedArray<T>(),
         m_binUtility(0)
        {}
        
        /**Constructor with std::vector and a BinUtility - reference counted, will delete objects at the end,
         if this deletion should be turned off, the boolean deletion should be switched to false
         the global position is given by object! */
        BinnedArray2D(const std::vector< std::pair< T, Vector3D > >& tclassvector, BinUtility* bingen) :
          BinnedArray<T>(),
          m_binUtility(bingen)
        {
            if (bingen){
                m_array = std::vector< std::vector< T > >(bingen->bins(1));
                for (size_t i=0; i < bingen->bins(1); ++i)
                    m_array[i] = std::vector< T >(bingen->bins(0),nullptr);
                // fill the Volume vector into the array
                size_t vecsize = tclassvector.size();
                for (size_t ivec = 0; ivec < vecsize; ++ivec) {
                    Vector3D currentGlobal(((tclassvector[ivec]).second));
                    if (bingen->inside(currentGlobal)){
                        // the object
                        T object = ((tclassvector)[ivec]).first;
                        // assign to the array store
                        ((m_array)[bingen->bin(currentGlobal,1)])[bingen->bin(currentGlobal,0)] = object;
                        // also fill the one-dimensional store - only if unique ones though
                        auto beginIter = BinnedArray<T>::m_arrayObjects.begin();
                        auto endIter   = BinnedArray<T>::m_arrayObjects.end();
                        if (std::find(beginIter,endIter,object) == endIter)
                            BinnedArray<T>::m_arrayObjects.push_back(object);
                    } else
                        throw "BinnedArray2D";
                }
            } else
              throw "BinnedArray2D";
        }
        
        /**Copy Constructor - copies only pointers !*/
        BinnedArray2D(const BinnedArray2D& barr) :
          BinnedArray<T>(barr),
          m_array(barr.m_array),
          m_binUtility(nullptr)
        {
            m_binUtility = (barr.m_binUtility) ? barr.m_binUtility->clone() : nullptr;
        }
        
        /**Assignment operator*/
        BinnedArray2D& operator=(const BinnedArray2D& barr)
        {
            if (this != &barr){
                delete m_binUtility;
                // now refill - copy
                m_binUtility = (barr.m_binUtility) ? barr.m_binUtility->clone() : 0;
                // array copying
                m_array = barr.m_array;
                BinnedArray<T>::m_arrayObjects = barr.m_arrayObjects;
            }
            return *this;
        }
        
        /** Implizit Constructor */
        BinnedArray2D* clone() const
        { return new BinnedArray2D(*this); }
        
        /**Virtual Destructor*/
        ~BinnedArray2D()
        {
            delete m_binUtility;
        }
        
        /** Returns the object in the array from a local position */
        T object(const Vector2D& lp) const final
        {
            return (((m_array)[m_binUtility->bin(lp, 1)]))[m_binUtility->bin(lp, 0)];
        }
        
        /** Returns the object in the array from a global position */
        T object(const Vector3D& gp) const final
        {
            return (((m_array)[m_binUtility->bin(gp, 1)]))[m_binUtility->bin(gp, 0)];
        }
        
        /** Returns the entry object in the array from a global position */
        T entryObject(const Vector3D& pos) const final
        {
            return (((m_array)[m_binUtility->entry(pos, 1)]))[m_binUtility->entry(pos, 0)];
        }
        
        /** Returns the objects in an odered fashion */
        const std::vector< std::vector< T > >& arrayObjectsOrdered() const { return m_array; }
        
        /** Return the BinUtility*/
        const BinUtility* binUtility() const { return(m_binUtility); }
        
    private:

        std::vector< std::vector< T > >   m_array;        //!< vector of pointers to the class T
        BinUtility*                       m_binUtility;   //!< binUtility for retrieving and filling the Array
        
    };
    
    
} // end of namespace Acts

#endif // ACTS_GEOMETRYUTILS_BINNEDARRAY2D_H

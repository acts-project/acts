///////////////////////////////////////////////////////////////////
// BinnedArray1D.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_BINNEDARRAY1D_H
#define ACTS_GEOMETRYUTILS_BINNEDARRAY1D_H 1

// Geometry module
#include "ACTS/Utilities/BinnedArray.h"
#include "ACTS/Utilities/BinUtility.h"
//STL
#include <vector>

namespace Acts {

/** @class BinnedArray1D

    1-dimensional binned array based on a sorting as given by the BinUtitlity.

    @author Andreas.Salzburger@cern.ch
   */

  template <class T> class BinnedArray1D : public BinnedArray<T> {

    public:
     /**Default Constructor - needed for inherited classes */
     BinnedArray1D():
       BinnedArray<T>(),
       m_binUtility(0)
     {}

     /** Constructor with std::vector< pair < object, position> > and a  BinUtility to sort the objects into the bin */
    BinnedArray1D(const std::vector< std::pair< T, Vector3D > >& tclassvector, BinUtility* bingen) throw (std::string) :
       BinnedArray<T>(),
       m_binUtility(bingen)
     {
        // prepare the binned Array
        if (bingen){
          // get the size of the object array and reserve accoridngly    
          size_t vecsize = tclassvector.size();
          BinnedArray<T>::m_arrayObjects = std::vector<T>(vecsize, nullptr);
          // fill into the new
          for (size_t ivec = 0; ivec < vecsize; ++ivec){
            Vector3D currentGlobal(((tclassvector[ivec]).second));
            if (bingen->inside(currentGlobal))
               BinnedArray<T>::m_arrayObjects[bingen->bin(currentGlobal,0)] = ((tclassvector)[ivec]).first;
            else 
               throw "Object outside bounds";  
          }
       }
     }

     /**Copy Constructor - copies only pointers !*/
    BinnedArray1D(const BinnedArray1D& barr) throw (std::string) :
       BinnedArray<T>(barr),
       m_binUtility(nullptr)
     {
        // clone the Bin utility for the copy constructor
         if (barr.m_binUtility)
             m_binUtility = barr.m_binUtility->clone();
        else 
           throw "No BinUtility provided";  
     }
      
     /**Assignment operator*/
     BinnedArray1D& operator=(const BinnedArray1D& barr)
     {
        if (this != &barr){
          // now refill
          m_binUtility = (barr.m_binUtility) ? barr.m_binUtility->clone() : nullptr;
          BinnedArray<T>::m_arrayObjects = barr.m_arrayObjects;
        }
        return *this;
     }
      
     /** Implicit Constructor */
     BinnedArray1D* clone() const final
     { return new BinnedArray1D(*this); }

     /**Virtual Destructor*/
     ~BinnedArray1D()
     {
        delete m_binUtility;
     }

     /** Returns the pointer to the templated class object from the BinnedArray,
         it returns nullptr if not defined */
     T object(const Vector2D& lp) const final
     { return (BinnedArray<T>::m_arrayObjects[m_binUtility->bin(lp, 0)]); }

     /** Returns the pointer to the templated class object from the BinnedArray
         it returns nullptr if not defined */
     T object(const Vector3D& gp) const final
     { return (BinnedArray<T>::m_arrayObjects[m_binUtility->bin(gp, 0)]); }

     /** Returns the pointer to the templated class object from the BinnedArray - entry point*/
     T entryObject(const Vector3D& gp) const final
     { return (BinnedArray<T>::m_arrayObjects[m_binUtility->entry(gp, 0)]); }

     /** Return the BinUtility*/
     const BinUtility* binUtility() const { return(m_binUtility); }

    private:
      BinUtility*    m_binUtility;   //!< binUtility for retrieving and filling the Array

  };

} // end of namespace Acts

#endif // ACTS_GEOMETRYUTILS_BINNEDARRAY1D_H

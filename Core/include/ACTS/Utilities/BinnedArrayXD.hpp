// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BinnedArray.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_UTILITIES_BINNEDARRAYXD_H
#define ACTS_UTILITIES_BINNEDARRAYXD_H 1

#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include <vector>
#include <iostream>

class MsgStream;

namespace Acts {

/// @class BinnedArray
///
/// Avoiding a map search, the templated BinnedArray class can help
/// ordereing geometrical objects by providing a dedicated BinUtility.
///
/// This can be 0D, 1D, 2D, and 3D
///
/// the type of Binning is given defined through the BinUtility   
template <class T> class BinnedArrayXD : public BinnedArray<T>
{

/// typedef the object and position for readability
typedef std::pair<T, Vector3D>  TAP;
  
public:
  /// Constructor for single object
  BinnedArrayXD(T object)
    : BinnedArray<T>()
    , m_objectGrid()
    , m_arrayObjects({object})
    , m_binUtility(nullptr)
  {
    /// fill the single object into the object grid
    m_objectGrid[0][0][0] = object;
  }

  /// Constructor with std::vector and a BinUtility 
  /// - fills the internal data structur
  BinnedArrayXD(const std::vector<TAP>& tapvector,
                std::unique_ptr<BinUtility> bu)
    : BinnedArray<T>()
    , m_objectGrid(bu->bins(2), std::vector< std::vector<T> >
                  (bu->bins(1), std::vector<T>(bu->bins(0), nullptr)))
    , m_arrayObjects()
    , m_binUtility(std::move(bu))
  {
    /// reserve the right amount of data
    m_arrayObjects.reserve(tapvector.size());
    /// loop over the object & position for ordering
    for (auto& tap : tapvector) {
      /// check for inside 
      if (m_binUtility->inside(tap.second)) {
        // butil to the array store - if the bingen
        // dimension is smaller 1,2 it will provide 0
        size_t bin0 = m_binUtility->bin(tap.second, 0);
        size_t bin1 = m_binUtility->bin(tap.second, 1);
        size_t bin2 = m_binUtility->bin(tap.second, 2);
        /// fill the data  
        m_objectGrid[bin2][bin1][bin0] = tap.first;
        /// fill the unique m_arrayObjects
        auto beginIter = m_arrayObjects.begin();
        auto endIter   = m_arrayObjects.end();
        if (std::find(beginIter, endIter, tap.first) == endIter)
            m_arrayObjects.push_back(tap.first); 
      } 
    }
  }

  /// Copy constructor
  /// - not allowed, use the same array
  BinnedArrayXD(const BinnedArrayXD<T>& barr) = delete;

  /// Assignment operator
  /// - not allowed, use the same array
  BinnedArrayXD&
  operator=(const BinnedArrayXD<T>& barr) = delete;

  /// Destructor
  ~BinnedArrayXD() {}
  
  /// Returns the object in the array from a local position
  ///
  /// @TODO check if we can change to triple return at once
  ///
  /// @param lposition is the local position for the bin search
  /// @return is the object in that bin
  T
  object(const Vector2D& lposition) const final
  {
    if (m_binUtility){
      size_t bdim = m_binUtility->dimensions();
      size_t bin2 = bdim > 2 ? m_binUtility->bin(lposition, 2)  :  0; 
      size_t bin1 = bdim > 1 ? m_binUtility->bin(lposition, 1)  :  0; 
      size_t bin0 =  m_binUtility->bin(lposition, 0);
      return m_objectGrid[bin2][bin1][bin0];      
    }
    return m_objectGrid[0][0][0];
  }

  /// Returns the object in the array from a global position 
  ///
  /// @TODO check if we can change to triple return at once
  ///
  /// @param position is the global position for the bin search
  /// @return is the object in that bin
  T
  object(const Vector3D& position) const final
  {
    if (m_binUtility){
      size_t bdim = m_binUtility->dimensions();
      size_t bin2 = bdim > 2 ? m_binUtility->bin(position, 2)  :  0; 
      size_t bin1 = bdim > 1 ? m_binUtility->bin(position, 1)  :  0; 
      size_t bin0 =  m_binUtility->bin(position, 0);
      return m_objectGrid[bin2][bin1][bin0];      
    }
    return m_objectGrid[0][0][0];
  }

  /// Return all unqiue object
  const std::vector<T>&
  arrayObjects() const final
  {
    return m_arrayObjects;
  };

  /// Return the BinUtility
  const BinUtility*
  binUtility() const final
  {
    return (m_binUtility.get());
  }

private:
  /// the data store - a 3D array at default
  std::vector< std::vector< std::vector < T > > > m_objectGrid; 
  /// Vector of unique Array objects
  std::vector<T>                                  m_arrayObjects; 
  /// binUtility for retrieving and filling the Array 
  std::unique_ptr<BinUtility>                     m_binUtility;   

};

}  // end of namespace Acts

#endif // ACTS_UTILITIES_BINNEDARRAYXD_H

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

#ifndef ACTS_UTILITIES_BINNEDARRAY_H
#define ACTS_UTILITIES_BINNEDARRAY_H 1

#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include <vector>

namespace Acts {

///  @class BinnedArray
/// 
/// Pure virtual base class for Binned Array to avoid map searches
/// - there is only one restriction:
///   T must be of pointer type in order to be initialized withh nullptr
///   and to allow for nullptr return type
/// 
/// - the BinnedArray is designed for 0D, 1D, 2D, and 3D binning  
template <class T>
class BinnedArray
{
public:
  /// Default Constructor - needed for inherited classes 
  BinnedArray() {}
  
  /// Virtual Destructor
  virtual ~BinnedArray() {}
  
  /// Returns the object in the associated bin according the local position  
  ///
  /// @param lposition is the local position for the object retrieval
  /// @return the object according to the estimated bin
  virtual T
  object(const Vector2D& lposition) const = 0;

  /// Returns the object in the associated bin according the local position  
  ///
  /// @param position is the global position for the object retrieval
  /// @return the object according to the estimated bin
  virtual T
  object(const Vector3D& position) const = 0;

  /// Return all unqiue object
  virtual const std::vector<T>&
  arrayObjects() const = 0;

  /// Return the BinUtility
  /// - if returned 0 it is a 0D array 
  virtual const BinUtility*
  binUtility() const = 0;

};

}  // end of namespace Acts

#endif  // ACTS_UTILITIES_BINNEDARRAY_H

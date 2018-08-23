// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BinnedArray.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <array>
#include <vector>
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Definitions.hpp"

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
  BinnedArray() = default;
  /// Virtual Destructor
  virtual ~BinnedArray() = default;
  /// Returns the object in the associated bin according the local position
  ///
  /// @param lposition is the local position for the object retrieval
  /// @param bins is the bin triple to filled
  ///
  /// @return the object according to the estimated bin
  virtual T
  object(const Vector2D& lposition, std::array<size_t, 3>& bins) const = 0;

  /// Same method without bins for backward compatibility
  ///
  /// @param lposition is the local position for finding the obect
  ///
  /// @return the object according to the estimated bin
  virtual T
  object(const Vector2D& lposition) const
  {
    std::array<size_t, 3> bins;
    return object(lposition, bins);
  }

  /// Returns the object in the associated bin according the local position
  ///
  /// @param position is the global position for the object retrieval
  /// @param bin is the bin triple filled
  ///
  /// @return the object according to the estimated bin
  virtual T
  object(const Vector3D& position, std::array<size_t, 3>& bin) const = 0;

  /// Same method without bins for backward compatibility
  ///
  /// @param position is the global position for the object finding
  ///
  /// @return the object according to the estimated bin
  virtual T
  object(const Vector3D& position) const
  {
    std::array<size_t, 3> bins;
    return object(position, bins);
  }

  /// Returns the object found through global position search
  /// and their neighbor objects
  ///
  /// @todo check if this needs connectivity directive
  ///
  /// @param bin is the binning
  ///
  /// @return a vector of unique objects
  virtual std::vector<T>
  objectCluster(const std::array<size_t, 3>& bin) const = 0;

  /// Return all unqiue object
  /// @note this is the accessor to the
  /// @return the vector of all array objects
  virtual const std::vector<T>&
  arrayObjects() const = 0;

  /// Return the object grid multiple entries are allowed
  /// @return the object grid
  virtual const std::vector<std::vector<std::vector<T>>>&
  objectGrid() const = 0;

  /// Return the BinUtility
  /// - if returned 0 it is a 0D array
  /// @return plain pointer to the bin utility
  virtual const BinUtility*
  binUtility() const = 0;
};

}  // namespace Acts

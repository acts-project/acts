// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BinnedArrayXD.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <array>
#include <iostream>
#include <vector>
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArray.hpp"

class MsgStream;

namespace Acts {

/// @class BinnedArrayXD
///
/// Avoiding a map search, the templated BinnedArray class can help
/// ordereing geometrical objects by providing a dedicated BinUtility.
///
/// This can be 0D, 1D, 2D, and 3D in regular binning
///
/// the type of Binning is given defined through the BinUtility
template <class T>
class BinnedArrayXD : public BinnedArray<T>
{
  /// typedef the object and position for readability
  using TAP = std::pair<T, Vector3D>;

public:
  /// Constructor for single object
  ///
  /// @tparam object is the single object
  BinnedArrayXD(T object)
    : BinnedArray<T>()
    , m_objectGrid(1,
                   std::vector<std::vector<T>>(1, std::vector<T>(1, nullptr)))
    , m_arrayObjects({object})
    , m_binUtility(nullptr)
  {
    /// fill the single object into the object grid
    m_objectGrid[0][0][0] = object;
  }

  /// Constructor with std::vector and a BinUtility
  /// - fills the internal data structur
  ///
  /// @param tapvector is a vector of object and binning position
  /// @param bu is the unique bin utility for this binned array
  BinnedArrayXD(const std::vector<TAP>&           tapvector,
                std::unique_ptr<const BinUtility> bu)
    : BinnedArray<T>()
    , m_objectGrid(
          bu->bins(2),
          std::vector<std::vector<T>>(bu->bins(1),
                                      std::vector<T>(bu->bins(0), nullptr)))
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
        auto bins = m_binUtility->binTriple(tap.second);
        /// fill the data
        m_objectGrid[bins[2]][bins[1]][bins[0]] = tap.first;
        /// fill the unique m_arrayObjects
        if (std::find(m_arrayObjects.begin(), m_arrayObjects.end(), tap.first)
            == m_arrayObjects.end()) {
          m_arrayObjects.push_back(tap.first);
        }
      }
    }
  }

  /// Constructor with a grid and a BinUtility
  ///
  /// @param grid is the prepared object grid
  /// @param bu is the unique bin utility for this binned array
  BinnedArrayXD(const std::vector<std::vector<std::vector<T>>>& grid,
                std::unique_ptr<const BinUtility>               bu)
    : BinnedArray<T>()
    , m_objectGrid(grid)
    , m_arrayObjects()
    , m_binUtility(std::move(bu))
  {
    // get the total dimension
    size_t objects
        = m_binUtility->bins(0) * m_binUtility->bins(1) * m_binUtility->bins(2);
    /// reserve the right amount of data
    m_arrayObjects.reserve(objects);
    /// loop over the object & position for ordering
    for (auto& o2 : m_objectGrid) {
      for (auto& o1 : o2) {
        for (auto& o0 : o1) {
          if (o0) {
            /// fill the unique m_arrayObjects
            if (std::find(m_arrayObjects.begin(), m_arrayObjects.end(), o0)
                == m_arrayObjects.end()) {
              m_arrayObjects.push_back(o0);
            }
          }
        }
      }
    }
  }

  /// Copy constructor
  /// - not allowed, use the same array
  BinnedArrayXD(const BinnedArrayXD<T>& barr) = delete;

  /// Assignment operator
  /// - not allowed, use the same array
  BinnedArrayXD&
  operator=(const BinnedArrayXD<T>& barr)
      = delete;

  /// Destructor
  ~BinnedArrayXD() override = default;
  /// Returns the object in the array from a local position
  ///
  /// @todo check if we can change to triple return at once
  ///
  /// @param lposition is the local position for the bin search
  /// @param bins is the bin triple filled during this access
  ///
  /// @return is the object in that bin
  T
  object(const Vector2D& lposition, std::array<size_t, 3>& bins) const final
  {
    if (m_binUtility) {
      size_t bdim = m_binUtility->dimensions();
      bins[2]     = bdim > 2 ? m_binUtility->bin(lposition, 2) : 0;
      bins[1]     = bdim > 1 ? m_binUtility->bin(lposition, 1) : 0;
      bins[0]     = m_binUtility->bin(lposition, 0);
      return m_objectGrid[bins[2]][bins[1]][bins[0]];
    }
    return m_objectGrid[0][0][0];
  }

  // satisfy overload / override
  T
  object(const Vector2D& lposition) const override
  {
    std::array<size_t, 3> bins;
    return object(lposition, bins);
  }

  /// Returns the object in the array from a global position
  ///
  /// @param position is the global position for the bin search
  /// @param bins is the bins triple filled during access
  ///
  /// @return is the object in that bin
  T
  object(const Vector3D& position, std::array<size_t, 3>& bins) const final
  {
    if (m_binUtility) {
      size_t bdim = m_binUtility->dimensions();
      bins[2]     = bdim > 2 ? m_binUtility->bin(position, 2) : 0;
      bins[1]     = bdim > 1 ? m_binUtility->bin(position, 1) : 0;
      bins[0]     = m_binUtility->bin(position, 0);
      return m_objectGrid[bins[2]][bins[1]][bins[0]];
    }
    return m_objectGrid[0][0][0];
  }

  // satisfy overload / override
  T
  object(const Vector3D& position) const override
  {
    std::array<size_t, 3> bins;
    return object(position, bins);
  }

  /// Return all unqiue object
  /// @return vector of unique array objects
  const std::vector<T>&
  arrayObjects() const final
  {
    return m_arrayObjects;
  }

  /// Return the object grid
  /// multiple entries are allowed and wanted
  /// @return internal object grid
  const std::vector<std::vector<std::vector<T>>>&
  objectGrid() const final
  {
    return m_objectGrid;
  }

  /// Returns the object according to the bin triple
  /// and their neighbour objects (if different)
  ///
  /// @param binTriple is the binning
  ///
  /// @return a vector of unique objects
  std::vector<T>
  objectCluster(const std::array<size_t, 3>& binTriple) const override
  {
    // prepare the return vector
    std::vector<T> rvector;
    // reference bin object to be excluded
    T bObject = m_objectGrid[binTriple[2]][binTriple[1]][binTriple[0]];
    // get the dimensions first
    size_t bdim = m_binUtility->dimensions();
    // avoiding code duplication
    std::vector<size_t> zerorange = {0};
    // 2D bin
    std::vector<size_t> bin2values = (bdim > 2)
        ? m_binUtility->binningData()[2].neighbourRange(binTriple[2])
        : zerorange;
    // 1D bin
    std::vector<size_t> bin1values = (bdim > 1)
        ? m_binUtility->binningData()[1].neighbourRange(binTriple[1])
        : zerorange;
    // 0D bin
    std::vector<size_t> bin0values
        = m_binUtility->binningData()[0].neighbourRange(binTriple[0]);

    // do the loop
    for (auto b2 : bin2values) {
      for (auto b1 : bin1values) {
        for (auto b0 : bin0values) {
          // get the object
          T object = m_objectGrid[b2][b1][b0];
          if (object && object != bObject
              && std::find(rvector.begin(), rvector.end(), object)
                  == rvector.end()) {
            rvector.push_back(object);
          }
        }
      }
    }
    // return the ones you found
    return rvector;
  }

  /// Return the BinUtility
  /// @return plain pointer to the bin utility of this array
  const BinUtility*
  binUtility() const final
  {
    return (m_binUtility.get());
  }

private:
  /// the data store - a 3D array at default
  std::vector<std::vector<std::vector<T>>> m_objectGrid;
  /// Vector of unique Array objects
  std::vector<T> m_arrayObjects;
  /// binUtility for retrieving and filling the Array
  std::unique_ptr<const BinUtility> m_binUtility;
};
}  // end of namespace Acts

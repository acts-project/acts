// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DefaultDetectorElementBase.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <memory>
#include <vector>
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class Surface;

/// @class DetectorElementBase
///
/// This is the default base class for all tracking detector elements
/// with read-out relevant information. It provides the minimal interface
/// for the Acts proxy mechanism for surfaces, i.e. surfaces in the
/// Tracking geometry representing actual detection devices
///
class DetectorElementBase
{
public:
  /// Constructor
  DetectorElementBase() = default;

  /// virtual Destructor
  virtual ~DetectorElementBase() = default;

  /// Return the transform for the Element proxy mechanism
  virtual const Transform3D&
  transform() const = 0;

  /// Return surface association
  virtual const Surface&
  surface() const = 0;

  /// Returns the thickness of the module
  /// @return double that indicates the thickness of the module
  virtual double
  thickness() const = 0;

  /// Fast access to bin members
  /// Bin members are elements that are in the same geometric binning cell,
  /// such, e.g. backside modules in a multiplet detector
  ///
  /// @return vector of DetectorElementBase pointers
  const std::vector<const DetectorElementBase*>&
  binmembers() const;

  /// Reigster the bin members
  /// Bin members are elements that are in the same geometric binning cell,
  /// such, e.g. backside modules in a doublet/triplet detector
  ///
  /// @param binmembers DetectorElementBase objects in the same cell
  void
  registerBinmembers(std::vector<const DetectorElementBase*>& binmembers);

  /// Fast access to neighbours
  /// Neighbours are elements that are in an neighbouring geometric binning
  /// cell, such, e.g. next in phi, next in eta modules
  ///
  /// @return vector of DetectorElementBase pointers
  const std::vector<const DetectorElementBase*>&
  neighbours() const;

  /// Reigster the neighbours
  /// Neighbours are elements that are in an neighbouring geometric binning
  /// cell, such, e.g. next in phi, next in eta modules
  ///
  /// @param neighbours are DetectorElementBase objects that are neighbours
  void
  registerNeighbours(std::vector<const DetectorElementBase*>& neighbours);

private:
  std::vector<const DetectorElementBase*> m_binmembers;
  std::vector<const DetectorElementBase*> m_neighbours;
};

inline const std::vector<const DetectorElementBase*>&
DetectorElementBase::binmembers() const
{
  return m_binmembers;
}

inline const std::vector<const DetectorElementBase*>&
DetectorElementBase::neighbours() const
{
  return m_neighbours;
}

inline void
DetectorElementBase::registerBinmembers(
    std::vector<const DetectorElementBase*>& binmembers)
{
  for (auto& bmember : binmembers) {
    // only fill if it's not yet registered
    if (find(m_binmembers.begin(), m_binmembers.end(), bmember)
        == m_binmembers.end()) {
      m_binmembers.push_back(bmember);
    }
  }
}

inline void
DetectorElementBase::registerNeighbours(
    std::vector<const DetectorElementBase*>& neighbours)
{
  for (auto& neighbour : neighbours) {
    // only fill if it's not yet registered
    if (find(m_neighbours.begin(), m_neighbours.end(), neighbour)
        == m_neighbours.end()) {
      m_neighbours.push_back(neighbour);
    }
  }
}

}  // end of namespace Acts

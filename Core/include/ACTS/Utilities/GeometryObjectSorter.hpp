// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GeometryObjectSorter.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYUTILS_GEOMETRYOBJECTSORTER_H
#define ACTS_GEOMETRYUTILS_GEOMETRYOBJECTSORTER_H 1

#include "Definitions.hpp"

namespace Acts {
/// @class GeometryObjectSorter
///
template <class T>
class GeometryObjectSorterT : public std::binary_function<T, T, bool>
{
public:
  /// Constructor from a binning value
  ///
  /// @param bValue is the value in which the binning is done
  /// @param transform is an optional transform to be performed
  GeometryObjectSorterT(BinningValue                 bValue,
                        std::shared_ptr<Transform3D> transform = nullptr)
    : m_binningValue(bValue), m_transform(transform)
  {
  }

  /// Comparison operator
  ///
  /// @tparam one first object
  /// @tparam two second object
  ///
  /// @return boolen indicator
  bool
  operator()(T one, T two) const
  {
    // get the pos one / pos two
    Vector3D posOne = m_transform
        ? m_transform->inverse() * one->binningPosition(m_binningValue)
        : one->binningPosition(m_binningValue);
    Vector3D posTwo = m_transform
        ? m_transform->inverse() * two->binningPosition(m_binningValue)
        : two->binningPosition(m_binningValue);

    // switch the binning value
    // - binX, binY, binZ, binR, binPhi, binRPhi, binH, binEta
    switch (m_binningValue) {
    // compare on x
    case binX: {
      return (posOne.x() < posTwo.x());
    }
    // compare on y
    case binY: {
      return (posOne.y() < posTwo.y());
    }
    // compare on z
    case binZ: {
      return (posOne.z() < posTwo.z());
    }
    // compare on r
    case binR: {
      return (posOne.perp() < posTwo.perp());
    }
    // compare on phi
    case binPhi: {
      return (posOne.phi() < posTwo.phi());
    }
    // compare on eta
    case binEta: {
      return (posOne.eta() < posTwo.eta());
    }
    // default for the moment
    default: {
      return (posOne.mag() < posTwo.mag());
    }
    }
  }

protected:
  BinningValue                 m_binningValue;
  std::shared_ptr<Transform3D> m_transform;
};
}

#endif  // ACTS_GEOMETRYUTILS_GEOMETRYOBJECTSORTER_H

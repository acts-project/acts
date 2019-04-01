// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ElementFraction.hpp, Acts project
///////////////////////////////////////////////////////////////////

/// We allow for a maximum of 256 elements (and no isotopes)
static const double s_oneOverUcharMax = 1. / double(UCHAR_MAX);

#pragma once
#include <algorithm>
#include <array>
#include <climits>
#include <iomanip>
#include <iostream>
#include <limits>

namespace Acts {

/// @class ElementFraction
///
/// This is a simple pair of the fractional component of an element
/// the first member is the element atomic charge to identify it uniquely
/// the second member is the fractional component
///
///
/// The elemnt fraction allows you to log with what fraction you have elements
/// in a merged material bin with many many bins. This is "A LOT" of
/// information, hence the data format has to be as small as possible
/// this serves only one purpose: multiple scattering and ionization loss does
/// not care about this at all, but nuclear interaction in the fast simulation
/// does.
///
/// uint8_t gives you 256 Elements (more than we need), with an accuracy of
/// 1./256, i.e. < 0.5 %. That's better than we can dream of parameterizing
/// hadronic interactions.
class ElementFraction
{
public:
  /// Default Constructor
  ElementFraction() = default;

  /// Constructor from arguments
  ///
  /// @param iz is the z value of the element as an unsigned int
  /// @param ifrac is the associated fraction of that element
  ElementFraction(unsigned int iz, float ifrac)
    : m_data{{(uint8_t)iz, (uint8_t)(ifrac * double(UCHAR_MAX))}}
  {
  }

  /// Constructor direct data
  ///
  /// @param data is the element fraction source object
  ElementFraction(const std::array<uint8_t, 2>& data) : m_data(data) {}

  /// Copy constructor
  ///
  /// @param ef is the element fraction source object
  ElementFraction(const ElementFraction& ef) = default;

  /// Copy move constructor
  ///
  /// @param ef is the element fraction source object
  ElementFraction(ElementFraction&& ef) = default;

  /// Assignment operator from base class
  ///
  /// @param ef is the element fraction source object
  ElementFraction&
  operator=(const ElementFraction& ef)
      = default;

  /// Assigment Move Operator
  ///
  /// @param ef is the element fraction source object
  ElementFraction&
  operator=(ElementFraction&& ef)
      = default;

  /// Access to the data itself
  const std::array<uint8_t, 2>&
  data() const
  {
    return m_data;
  }

  /// Return in a nice format
  /// @return casts back to an unsigned integer
  unsigned int
  element() const
  {
    return static_cast<unsigned int>(m_data[0]);
  }

  /// Return in a nice format
  /// @return casts char to an unsigned int and then into double
  double
  fraction() const
  {
    return (static_cast<unsigned int>(m_data[1]) * s_oneOverUcharMax);
  }

  /// Define the equality operator
  /// @param ef is the source ElementFraction for comparison
  bool
  operator==(const ElementFraction& ef) const
  {
    return (m_data[0] == ef.m_data[0] && m_data[1] == ef.m_data[1]);
  }

  /// Define smaller operator for sorting
  /// we always sort by fraction for fastest access to the
  /// most probable fraction
  ///
  /// @param ef is the source ElementFraction for comparison
  bool
  operator<(const ElementFraction& ef) const
  {
    return (m_data[1] < ef.m_data[1]);
  }

private:
  std::array<uint8_t, 2> m_data = {{0, 0}};  //!< the data component
};
}
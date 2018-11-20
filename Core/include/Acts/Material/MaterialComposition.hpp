// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialComposition.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "Acts/Material/ElementFraction.hpp"

namespace Acts {

/// @class MaterialComposition
///
/// This helper struct allows to create a material composition
/// as terms of element fraction objects
class MaterialComposition
{
public:
  /// Default constructor
  MaterialComposition() = default;

  /// Destructor
  ~MaterialComposition() = default;

  /// Constructor from vector of pairs
  ///
  /// It does perform a rescaling sucht that the fractions add up
  /// to one within tolerance
  ///
  /// @param efracs are the element fractions
  MaterialComposition(const std::vector<ElementFraction>& efracs)
  {
    m_elements = efracs;
    std::sort(m_elements.begin(), m_elements.end());
  }

  /// Copy Constructor
  ///
  /// @param mc is the element fraction vector
  MaterialComposition(const MaterialComposition& mc) = default;

  /// Copy Move constructor
  ///
  /// @param mc is the element fraction vector
  MaterialComposition(MaterialComposition&& mc) = default;

  /// Assignment operator
  ///
  /// @param mc is the source object
  MaterialComposition&
  operator=(const MaterialComposition& mc)
      = default;

  /// Assignment move operator
  ///
  /// @param mc is the source object
  MaterialComposition&
  operator=(MaterialComposition&& mc)
      = default;

  /// Access to the elements themselves
  const std::vector<ElementFraction>&
  elements() const
  {
    return m_elements;
  }

  /// Boolean operator to indicate if this is empty
  operator bool() const { return !empty(); }

  /// How many elements you have
  size_t
  size() const
  {
    return m_elements.size();
  }

  /// Check if empty
  bool
  empty() const
  {
    return m_elements.empty();
  }

  /// Euality operator
  bool
  operator==(const MaterialComposition& mc) const
  {
    return (mc.m_elements == m_elements);
  }

private:
  std::vector<ElementFraction> m_elements = {};
};
}

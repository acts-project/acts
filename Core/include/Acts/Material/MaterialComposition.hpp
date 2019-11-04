// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

namespace Acts {

/// Memory-efficient storage of the relative fraction of an element.
///
/// This can be used to define materials that are compounds of multiple elements
/// with varying fractions. The element is identified by its atomic number
/// stored as a single byte (allows up to 256 elements; more than we need).
/// Its fraction is also stored as a a single byte with values between 0 and
/// 255. This gives an accuracy of 1/256 ~ 0.5 %.
///
/// The element fraction allows you to store element composition in merged
/// material in a merged materials with a large number of bins. Depending on the
/// detector and the description granularity this can be a lot of information
/// and thus requires the reduced memory footprint. This is really only needed
/// for nuclear interaction in the fast simulation where the reduced fractional
/// accuracy is not a problem. It should be much better than the parametrization
/// uncertainty for hadronic interactions.
class ElementFraction {
 public:
  /// Construct from atomic number and relative fraction.
  ///
  /// @param e is the atomic number of the element
  /// @param f is the relative fraction and must be a value in [0,1]
  constexpr ElementFraction(unsigned int e, float f)
      : m_element(static_cast<uint8_t>(e)),
        m_fraction(static_cast<uint8_t>(f * UINT8_MAX)) {
    assert((0u < e) and ("The atomic number must be positive"));
    assert((0.0f <= f) and (f <= 1.0f) and
           "Relative fraction must be in [0,1]");
  }
  /// Construct from atomic number and integer weight.
  ///
  /// @param e is the atomic number of the element
  /// @param w is the integer weight and must be a value in [0,256)
  constexpr explicit ElementFraction(unsigned int e, unsigned int w)
      : m_element(static_cast<uint8_t>(e)),
        m_fraction(static_cast<uint8_t>(w)) {
    assert((0u < e) and ("The atomic number must be positive"));
    assert((w < 256u) and "Integer weight must be in [0,256)");
  }

  /// Must always be created with valid data.
  ElementFraction() = delete;
  ElementFraction(ElementFraction&&) = default;
  ElementFraction(const ElementFraction&) = default;
  ~ElementFraction() = default;
  ElementFraction& operator=(ElementFraction&&) = default;
  ElementFraction& operator=(const ElementFraction&) = default;

  /// The element atomic number.
  constexpr uint8_t element() const { return m_element; }
  /// The relative fraction of this element.
  constexpr float fraction() const {
    return static_cast<float>(m_fraction) / UINT8_MAX;
  }

 private:
  // element atomic number
  uint8_t m_element;
  // element fraction in the compound scaled to the [0,256) range.
  uint8_t m_fraction;

  friend constexpr bool operator==(ElementFraction lhs, ElementFraction rhs) {
    return (lhs.m_fraction == rhs.m_fraction) and
           (lhs.m_element == rhs.m_element);
  }
  /// Sort by fraction for fastest access to the most probable element.
  friend constexpr bool operator<(ElementFraction lhs, ElementFraction rhs) {
    return lhs.m_fraction < rhs.m_fraction;
  }
};

/// @class MaterialComposition
///
/// This helper struct allows to create a material composition
/// as terms of element fraction objects
class MaterialComposition {
 public:
  MaterialComposition() = default;

  /// Destructor
  ~MaterialComposition() = default;

  /// Constructor from vector of pairs
  ///
  /// It does perform a rescaling sucht that the fractions add up
  /// to one within tolerance
  ///
  /// @param efracs are the element fractions
  MaterialComposition(const std::vector<ElementFraction>& efracs) {
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
  MaterialComposition& operator=(const MaterialComposition& mc) = default;

  /// Assignment move operator
  ///
  /// @param mc is the source object
  MaterialComposition& operator=(MaterialComposition&& mc) = default;

  /// Access to the elements themselves
  const std::vector<ElementFraction>& elements() const { return m_elements; }

  /// Boolean operator to indicate if this is empty
  operator bool() const { return !empty(); }

  /// How many elements you have
  size_t size() const { return m_elements.size(); }

  /// Check if empty
  bool empty() const { return m_elements.empty(); }

  /// Euality operator
  bool operator==(const MaterialComposition& mc) const {
    return (mc.m_elements == m_elements);
  }

 private:
  std::vector<ElementFraction> m_elements = {};
};

}  // namespace Acts

// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iosfwd>
#include <limits>

#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// Material description for interactions with matter.
///
/// The following parameters are used to specify the material and its
/// interactions with traversing particles:
///
/// - radiation length X0 (native length units)
/// - nuclear interaction length L0 (native length units)
/// - relative atomic mass Ar (unitless number)
/// - nuclear charge number Z (elementary charge e)
/// - mass density rho (native mass unit / (native length unit)³)
///
/// The parameters can be effective or average parameters when e.g. a mixture
/// of materials is described.
class Material {
 public:
  /// Construct a vacuum representation.
  Material() = default;
  /// Construct from material parameters.
  ///
  /// @param X0_  is the radiation length
  /// @param L0_  is the nuclear interaction length
  /// @param Ar_  is the relative atomic mass
  /// @param Z_   is the atomic number
  /// @param rho_ is the mass density
  constexpr Material(float X0_, float L0_, float Ar_, float Z_, float rho_)
      : m_x0(X0_), m_l0(L0_), m_ar(Ar_), m_z(Z_), m_rho(rho_) {}
  /// Construct from an encoded parameters vector.
  Material(const ActsVectorF<5>& parameters);

  Material(Material&& mat) = default;
  Material(const Material& mat) = default;
  ~Material() = default;
  Material& operator=(Material&& mat) = default;
  Material& operator=(const Material& mat) = default;

  /// Check if the material is valid, i.e. it is not vacuum.
  constexpr operator bool() const { return 0.0f < m_ar; }

  /// Return the radition length. Infinity in case of vacuum.
  constexpr float X0() const { return m_x0; }
  /// Return the nuclear interaction length. Infinity in case of vacuum.
  constexpr float L0() const { return m_l0; }
  /// Return the relative atomic mass.
  constexpr float Ar() const { return m_ar; }
  /// Return the nuclear charge number.
  constexpr float Z() const { return m_z; }
  /// Return the mass density.
  constexpr float massDensity() const { return m_rho; }
  /// Return the molar electron density in mol / (native length unit)³.
  ///
  /// Use mol instead of the real number of electrons to avoid large numbers
  /// which could result in numerical instabilities somewhere else.
  float molarElectronDensity() const;
  /// Return the mean electron excitation energy.
  float meanExcitationEnergy() const;

  /// Encode the properties into an opaque parameters vector.
  ActsVectorF<5> classificationNumbers() const;

 private:
  float m_x0 = std::numeric_limits<float>::infinity();
  float m_l0 = std::numeric_limits<float>::infinity();
  float m_ar = 0.0f;
  float m_z = 0.0f;
  float m_rho = 0.0f;

  friend constexpr bool operator==(const Material& lhs, const Material& rhs) {
    return (lhs.m_x0 == rhs.m_x0) and (lhs.m_l0 == rhs.m_l0) and
           (lhs.m_ar == rhs.m_ar) and (lhs.m_z == rhs.m_z) and
           (lhs.m_rho == rhs.m_rho);
  }
  friend constexpr bool operator!=(const Material& lhs, const Material& rhs) {
    return !(lhs == rhs);
  }
};

std::ostream& operator<<(std::ostream& os, const Material& material);

}  // namespace Acts

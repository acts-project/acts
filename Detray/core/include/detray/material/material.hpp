// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// Detray include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/material/detail/density_effect_data.hpp"
#include "detray/utils/hash_combine.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <cstddef>
#include <functional>
#include <ratio>
#include <sstream>

namespace detray {

/// Material composition properties and density data used by interactions.
template <concepts::scalar scalar_t, typename R = std::ratio<1, 1>>
struct material {
  /// Compile-time fraction represented by this material component.
  using ratio = R;
  /// Scalar type for material properties.
  using scalar_type = scalar_t;

  constexpr material() = default;

  /// Construct from material parameters
  /// @param x0 radiation length
  /// @param l0 nuclear interaction length
  /// @param ar relative atomic mass
  /// @param z nuclear charge number
  /// @param mass_rho mass density
  DETRAY_HOST_DEVICE
  constexpr material(const scalar_type x0, const scalar_type l0,
                     const scalar_type ar, const scalar_type z,
                     const scalar_type mass_rho)
      : m_x0(x0), m_l0(l0), m_ar(ar), m_z(z), m_mass_rho(mass_rho) {
    m_molar_rho = mass_to_molar_density(ar, mass_rho);
  }

  /// Construct from material parameters with density effect data
  /// @param x0 radiation length
  /// @param l0 nuclear interaction length
  /// @param ar relative atomic mass
  /// @param z nuclear charge number
  /// @param mass_rho mass density
  /// @param d_a density-effect fit parameter a
  /// @param d_m density-effect fit parameter m
  /// @param d_X0 lower density-effect threshold
  /// @param d_X1 upper density-effect threshold
  /// @param d_I mean excitation energy
  /// @param d_nC density-effect normalization
  /// @param d_delta0 low-energy density-effect correction
  DETRAY_HOST_DEVICE
  constexpr material(const scalar_type x0, const scalar_type l0,
                     const scalar_type ar, const scalar_type z,
                     const scalar_type mass_rho, const scalar_type d_a,
                     const scalar_type d_m, const scalar_type d_X0,
                     const scalar_type d_X1, const scalar_type d_I,
                     const scalar_type d_nC, const scalar_type d_delta0)
      : m_x0(x0),
        m_l0(l0),
        m_ar(ar),
        m_z(z),
        m_mass_rho(mass_rho),
        m_density(d_a, d_m, d_X0, d_X1, d_I, d_nC, d_delta0),
        m_has_density_effect_data(true) {
    m_molar_rho = mass_to_molar_density(ar, mass_rho);
  }

  /// Equality operator
  ///
  /// @param rhs is the right hand side to be compared to
  /// @return whether radiation length, interaction length, atomic mass, and charge agree
  DETRAY_HOST_DEVICE
  constexpr bool operator==(const material<scalar_type> &rhs) const {
    // Don't check dependent parameters
    assert(m_molar_rho == mass_to_molar_density(m_ar, m_mass_rho));
    assert(rhs.molar_density() ==
           mass_to_molar_density(rhs.Ar(), rhs.mass_density()));

    return (m_x0 == rhs.X0() && m_l0 == rhs.L0() && m_ar == rhs.Ar() &&
            m_z == rhs.Z() && m_mass_rho == rhs.mass_density() &&
            m_has_density_effect_data == rhs.has_density_effect_data() &&
            m_density == rhs.density_effect_data());
  }

  /// Equality operator
  ///
  /// @param rhs is the right hand side to be compared to
  /// @return whether radiation length, interaction length, atomic mass, and charge agree
  DETRAY_HOST_DEVICE
  constexpr bool operator!=(const material<scalar_type> &rhs) const {
    return !(*this == rhs);
  }

  /// Boolean operator
  /// @returns whether the material contains valid values
  DETRAY_HOST_DEVICE
  constexpr bool is_valid() const {
    if (m_x0 == detail::invalid_value<scalar_type>()) {
      return false;
    }
    // Cannot happen independently from X0 not being invalid
    assert(m_l0 != detail::invalid_value<scalar_type>());
    assert(m_ar != 0.f);
    assert(m_z != 0.f);
    assert(m_mass_rho != 0.f);
    assert(m_molar_rho);

    return true;
  }

  /// @returns the radition length. Infinity in case of vacuum.
  DETRAY_HOST_DEVICE
  constexpr scalar_type X0() const { return m_x0; }
  /// @returns the nuclear interaction length. Infinity in case of vacuum.
  DETRAY_HOST_DEVICE
  constexpr scalar_type L0() const { return m_l0; }
  /// @returns the relative atomic mass.
  DETRAY_HOST_DEVICE
  constexpr scalar_type Ar() const { return m_ar; }
  /// @returns the nuclear charge number.
  DETRAY_HOST_DEVICE
  constexpr scalar_type Z() const { return m_z; }
  /// @returns the mass density.
  DETRAY_HOST_DEVICE
  constexpr scalar_type mass_density() const { return m_mass_rho; }
  /// @returns the molar density.
  /// @brief mass_density / A
  DETRAY_HOST_DEVICE
  constexpr scalar_type molar_density() const { return m_molar_rho; }
  /// @returns the molar electron density.
  /// @brief (Z/A) * mass_density
  DETRAY_HOST_DEVICE
  constexpr scalar_type molar_electron_density() const {
    return m_z * m_molar_rho;
  }

  /// @returns the density effect data
  DETRAY_HOST_DEVICE
  constexpr const detail::density_effect_data<scalar_type> &
  density_effect_data() const {
    return m_density;
  }

  /// @returns the (Approximated) mean excitation energy
  DETRAY_HOST_DEVICE
  scalar_type mean_excitation_energy() const {
    if (!m_has_density_effect_data) {
      // use approximative computation as defined in ATL-SOFT-PUB-2008-003
      return 16.f * unit<scalar_type>::eV *
             math::pow(m_z, static_cast<scalar_type>(0.9));
    } else {
      return m_density.get_mean_excitation_energy();
    }
  }

  /// Return this component's compile-time mass fraction.
  /// @return ratio numerator divided by denominator, or an invalid value for a zero denominator
  DETRAY_HOST_DEVICE
  constexpr scalar_type fraction() const {
    if constexpr (ratio::num == 0) {
      return 0.f;
    } else if constexpr (ratio::den == 0) {
      return detail::invalid_value<scalar_type>();
    } else {
      return static_cast<scalar_type>(ratio::num) /
             static_cast<scalar_type>(ratio::den);
    }
  }

  /// @returns a string that contains the material details
  DETRAY_HOST
  std::string to_string() const {
    std::stringstream strm;
    if (m_ar <= 0.f) {
      strm << "vacuum";
      return strm.str();
    }
    strm << "material: ";
    strm << " X0 = " << m_x0 << " mm";
    strm << " | L0 = " << m_l0 << " mm";
    strm << " | Z = " << m_z;
    strm << " | rho = " << m_mass_rho << " GeV/mm^3";

    return strm.str();
  }

  /// @returns a string stream that prints the material details
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &os, const material &mat) {
    os << mat.to_string();
    return os;
  }

  /// @returns true if density effect data is present
  DETRAY_HOST_DEVICE
  bool has_density_effect_data() const { return m_has_density_effect_data; }

 protected:
  /// Set the radition length.
  /// @param x0 new radiation length
  DETRAY_HOST_DEVICE
  constexpr void set_X0(scalar_type x0) { m_x0 = x0; }
  /// Set the nuclear interaction length.
  /// @param l0 new nuclear interaction length
  DETRAY_HOST_DEVICE
  constexpr void set_L0(scalar_type l0) { m_l0 = l0; }
  /// Set the relative atomic mass.
  /// @param ar new relative atomic mass
  DETRAY_HOST_DEVICE
  constexpr void set_Ar(scalar_type ar) { m_ar = ar; }
  /// Set the nuclear charge number.
  /// @param z new nuclear charge number
  DETRAY_HOST_DEVICE
  constexpr void set_Z(scalar_type z) { m_z = z; }
  /// Set the mass density.
  /// @param mass_rho new mass density
  DETRAY_HOST_DEVICE
  constexpr void set_mass_density(scalar_type mass_rho) {
    m_mass_rho = mass_rho;
  }
  /// Set the molar density.
  /// @param molar_rho new molar density
  DETRAY_HOST_DEVICE
  constexpr void set_molar_density(scalar_type molar_rho) {
    m_molar_rho = molar_rho;
  }
  /// Convert mass density to molar density using relative atomic mass.
  /// @param ar relative atomic mass
  /// @param mass_rho mass density
  /// @return molar density
  DETRAY_HOST_DEVICE
  constexpr scalar_type mass_to_molar_density(double ar,
                                              double mass_rho) const {
    if (mass_rho == 0.) {
      return 0.f;
    }

    const double molar_mass{ar * unit<double>::g / unit<double>::mol};

    return static_cast<scalar_type>(mass_rho / molar_mass);
  }

 private:
  // Material properties
  scalar_type m_x0 = detail::invalid_value<scalar_type>();
  scalar_type m_l0 = detail::invalid_value<scalar_type>();
  scalar_type m_ar = 0.f;
  scalar_type m_z = 0.f;
  scalar_type m_mass_rho = 0.f;
  scalar_type m_molar_rho = 0.f;
  detail::density_effect_data<scalar_type> m_density = {};
  bool m_has_density_effect_data = false;
};

/// Macro for declaring the predefined materials (w/o Density effect data)
/// @param MATERIAL_NAME generated material type name
/// @param X0 radiation length
/// @param L0 nuclear interaction length
/// @param Ar relative atomic mass
/// @param Z nuclear charge number
/// @param Rho mass density
#define DETRAY_DECLARE_MATERIAL(MATERIAL_NAME, X0, L0, Ar, Z, Rho)    \
  template <concepts::scalar scalar_t, typename R = std::ratio<1, 1>> \
  struct MATERIAL_NAME final : public material<scalar_t, R> {         \
    /** Material storage and property accessors. */                   \
    using base_type = material<scalar_t, R>;                          \
    using base_type::base_type;                                       \
    DETRAY_HOST_DEVICE                                                \
    constexpr MATERIAL_NAME() : base_type(X0, L0, Ar, Z, Rho) {}      \
  }

/// !EXPERIMENTAL!
/// Macro for declaring the predefined materials (with Density effect data)
/// @param MATERIAL_NAME generated material type name
/// @param X0 radiation length
/// @param L0 nuclear interaction length
/// @param Ar relative atomic mass
/// @param Z nuclear charge number
/// @param Rho mass density
/// @param Density0 density-effect fit parameter a
/// @param Density1 density-effect fit parameter m
/// @param Density2 lower density-effect threshold
/// @param Density3 upper density-effect threshold
/// @param Density4 mean excitation energy
/// @param Density5 density-effect normalization
/// @param Density6 low-energy density-effect correction
#define DETRAY_DECLARE_MATERIAL_WITH_DED(                                      \
    MATERIAL_NAME, X0, L0, Ar, Z, Rho, Density0, Density1, Density2, Density3, \
    Density4, Density5, Density6)                                              \
  template <concepts::scalar scalar_t, typename R = std::ratio<1, 1>>          \
  struct MATERIAL_NAME final : public material<scalar_t, R> {                  \
    /** Material storage and property accessors. */                            \
    using base_type = material<scalar_t, R>;                                   \
    using base_type::base_type;                                                \
    DETRAY_HOST_DEVICE                                                         \
    constexpr MATERIAL_NAME()                                                  \
        : base_type(X0, L0, Ar, Z, Rho, Density0, Density1, Density2,          \
                    Density3, Density4, Density5, Density6) {}                 \
  }

}  // namespace detray

namespace std {

/// Specialize std::hash, so that the material can be used e.g. in an
/// unordered_map or for the deduplication of the detector material
///
/// Uses the independent material parameters, which are compared by the
/// equality operator of the material: Identical materials have the same hash.
///
/// @tparam scalar_t the scalar type of the material parameters
/// @tparam R the ratio of the material in a mixture
template <detray::concepts::scalar scalar_t, typename R>
struct hash<detray::material<scalar_t, R>> {
  /// @param mat the material
  /// @returns the hash value of @p mat
  std::size_t operator()(
      const detray::material<scalar_t, R> &mat) const noexcept {
    std::size_t seed{0u};
    seed = detray::hash_combine(seed, mat.X0());
    seed = detray::hash_combine(seed, mat.L0());
    seed = detray::hash_combine(seed, mat.Ar());
    seed = detray::hash_combine(seed, mat.Z());
    seed = detray::hash_combine(seed, mat.mass_density());
    return seed;
  }
};

}  // namespace std

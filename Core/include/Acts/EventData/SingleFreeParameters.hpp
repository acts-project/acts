// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iomanip>
#include <ostream>

#include "Acts/EventData/detail/coordinate_transformations.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class SingleFreeParameters
///
/// @brief Container class for free parameters
///
/// This is a base class for neutral and charged free parameters. All parameters
/// and the corresponding covariance matrix is stored in global coordinates. It
/// is assumed that the order of the parameters (and thereby the order of the
/// entries of the covariance as well) is given as (position_x, position_y,
/// position_z, time, direction_x, direction_y, direction_z, charge /
/// |momentum|).
/// @tparam ChargePolicy Parameter that describes if the particle is charged or
/// neutral
/// @note It is assumed that a charged particle has a charge of +/-1
template <class ChargePolicy>
class SingleFreeParameters {
  static_assert(std::is_same<ChargePolicy, ChargedPolicy>::value or
                    std::is_same<ChargePolicy, NeutralPolicy>::value,
                "ChargePolicy must either be 'Acts::ChargedPolicy' or "
                "'Acts::NeutralPolicy");

 public:
  /// Public typedefs
  /// Type of covariance matrix
  using CovMatrix_t = FreeSymMatrix;

  /// @brief Default virtual destructor
  ~SingleFreeParameters() = default;

  /// @brief Standard constructor for track parameters of charged particles
  ///
  /// @tparam T Type of the charge policy (ChargedPolicy)
  /// @param [in] cov The covariance matrix
  /// @param [in] parValues Vector with parameter values
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  SingleFreeParameters(std::optional<CovMatrix_t> cov,
                       const FreeVector& parValues)
      : m_parameters(parValues),
        m_oChargePolicy((0 < parValues(eFreeQOverP)) ? 1. : -1.),
        m_covariance(std::move(cov)) {}

  /// @brief Standard constructor for track parameters of neutral particles
  ///
  /// @tparam T Type of the charge policy (NeutralPolicy)
  /// @param [in] cov The covariance matrix
  /// @param [in] parValues Vector with parameter values
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  SingleFreeParameters(std::optional<CovMatrix_t> cov,
                       const FreeVector& parValues)
      : m_parameters(parValues),
        m_oChargePolicy(),
        m_covariance(std::move(cov)) {}

  /// @brief Copy assignment operator
  ///
  /// @param [in] rhs Object to be copied
  ///
  /// @return The assigned-to object `*this`
  SingleFreeParameters<ChargePolicy>& operator=(
      const SingleFreeParameters<ChargePolicy>& rhs) {
    // Check for self-assignment
    if (this != &rhs) {
      m_oChargePolicy = rhs.m_oChargePolicy;
      m_parameters = rhs.m_parameters;
      m_covariance = rhs.m_covariance;
    }
    return *this;
  }

  /// @brief Move assignment operator
  ///
  /// @param [in] rhs object to be movied into `*this`
  ///
  /// @return The assigned-to object `*this`
  SingleFreeParameters<ChargePolicy>& operator=(
      SingleFreeParameters<ChargePolicy>&& rhs) {
    // Check for self-assignment
    if (this != &rhs) {
      m_oChargePolicy = std::move(rhs.m_oChargePolicy);
      m_parameters = std::move(rhs.m_parameters);
      m_covariance = std::move(rhs.m_covariance);
    }
    return *this;
  }

  /// @brief Default copy constructor
  ///
  /// @param [in] copy The object to copy from
  SingleFreeParameters(const SingleFreeParameters<ChargePolicy>& copy)
      : m_parameters(copy.m_parameters),
        m_oChargePolicy(copy.m_oChargePolicy),
        m_covariance(copy.m_covariance) {}

  /// @brief Default move constructor
  ///
  /// @param [in] copy The object to move from
  SingleFreeParameters(SingleFreeParameters<ChargePolicy>&& copy) {
    this->operator=(
        std::forward<const SingleFreeParameters<ChargePolicy>>(copy));
  }

  /// @brief Access all parameters
  ///
  /// @return Vector containing the store parameters
  FreeVector parameters() const { return m_parameters; }

  /// @brief Access to a single parameter
  ///
  /// @tparam par Identifier of the parameter index which will be retrieved
  ///
  /// @return Value of the requested parameter
  template <unsigned int par,
            std::enable_if_t<par<eFreeParametersSize, int> = 0> ParValue_t get()
                const {
    return m_parameters(par);
  }

  /// @brief Access track parameter uncertainty
  ///
  /// @tparam par Identifier of the parameter uncertainty index which will be
  /// retrieved
  ///
  /// @return Value of the requested parameter uncertainty
  template <unsigned int par,
            std::enable_if_t<par<eFreeParametersSize, int> = 0> ParValue_t
                uncertainty() const {
    return std::sqrt(m_covariance->coeff(par, par));
  }

  /// @brief Access covariance matrix of track parameters
  ///
  /// @note The ownership of the covariance matrix is @b not transferred with
  /// this call.
  ///
  /// @return Raw pointer to covariance matrix (can be a nullptr)
  ///
  /// @sa ParameterSet::getCovariance
  const std::optional<CovMatrix_t>& covariance() const { return m_covariance; }

  /// @brief access position in global coordinate system
  ///
  /// @return 3D vector with global position
  Vector3D position() const { return m_parameters.template head<3>(); }

  /// @brief access momentum in global coordinate system
  ///
  /// @return 3D vector with global momentum
  Vector3D momentum() const {
    return m_parameters.template segment<3>(4) / std::abs(get<7>());
  }

  /// @brief retrieve electric charge
  ///
  /// @return value of electric charge
  double charge() const { return m_oChargePolicy.getCharge(); }

  /// @brief retrieve time
  ///
  /// @return value of time
  double time() const { return m_parameters(3); }

  /// @brief Equality operator
  ///
  /// @param [in] rhs Object to compare `*this` to
  ///
  /// @return Boolean value whether the objects can be casted into each other
  /// and the content of the member variables is the same
  bool operator==(const SingleFreeParameters& rhs) const {
    auto casted = dynamic_cast<decltype(this)>(&rhs);
    if (!casted) {
      return false;
    }

    // Both have covariance matrices set
    if ((m_covariance.has_value() && casted->m_covariance.has_value()) &&
        (*m_covariance != *casted->m_covariance)) {
      return false;
    }
    // Only one has a covariance matrix set
    if ((m_covariance.has_value() && !casted->m_covariance.has_value()) ||
        (!m_covariance.has_value() && casted->m_covariance.has_value())) {
      return false;
    }

    return (m_oChargePolicy == casted->m_oChargePolicy &&
            m_parameters == casted->m_parameters);
  }

  /// @brief inequality operator
  ///
  /// @return `not (*this == rhs)`
  bool operator!=(const SingleFreeParameters& rhs) const {
    return !(*this == rhs);
  }

  /// @brief Update of the parameterisation
  ///
  /// @tparam par The parameter index
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param newValue The new updaed value
  ///
  /// @note The context is not used here but makes the API consistent with @c
  /// SingleCurvilinearTrackParameters and @c SingleBoundTrackParameters
  template <unsigned int par,
            std::enable_if_t<par<eFreeParametersSize, int> = 0> void set(
                const GeometryContext& /*gctx*/, ParValue_t newValue) {
    m_parameters(par) = newValue;
  }

  /// @brief Print information to output stream
  ///
  /// @param [in, out] sl The output stream
  ///
  /// @return The modified output stream object @p sl
  std::ostream& print(std::ostream& sl) const {
    // Set stream output format
    auto old_precision = sl.precision(7);
    auto old_flags = sl.setf(std::ios::fixed);

    // Fill stream with content
    sl << " * FreeTrackParameters: ";
    sl << parameters().transpose() << std::endl;
    sl << " * charge: " << charge() << std::endl;
    if (covariance().has_value()) {
      sl << " * covariance matrix:\n" << *covariance() << std::endl;
    } else {
      sl << " * no covariance matrix stored" << std::endl;
    }

    // Reset stream format
    sl.precision(old_precision);
    sl.setf(old_flags);

    return sl;
  }

  /// @brief Output stream operator
  ///
  /// Prints information about this object to the output stream
  /// @param [in, out] out The output stream
  /// @param [in] sfp The object that will be printed
  ///
  /// @return Modified output stream object
  friend std::ostream& operator<<(std::ostream& out,
                                  const SingleFreeParameters& sfp) {
    sfp.print(out);
    return out;
  }

 private:
  FreeVector m_parameters;       ///< Parameter vector
  ChargePolicy m_oChargePolicy;  ///< charge policy object distinguishing
                                 /// between charged and neutral tracks
  std::optional<CovMatrix_t> m_covariance;  ///< Covariance matrix
};

}  // namespace Acts

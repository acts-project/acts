// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iomanip>
#include <ostream>

#include "Acts/EventData/ParameterSet.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// Free track parameters not bound to a surface.
///
/// This is a base class for neutral and charged free parameters. All parameters
/// and the corresponding covariance matrix is stored in global coordinates. It
/// is assumed that the order of the parameters (and thereby the order of the
/// entries of the covariance as well) is given as (position_x, position_y,
/// position_z, time, direction_x, direction_y, direction_z, charge /
/// |momentum|).
///
/// @tparam ChargePolicy Selection type if the particle is charged or neutral
/// @note It is assumed that a charged particle has a charge of +/-1
template <class ChargePolicy>
class SingleFreeTrackParameters {
  static_assert(std::is_same<ChargePolicy, ChargedPolicy>::value or
                    std::is_same<ChargePolicy, NeutralPolicy>::value,
                "ChargePolicy must either be 'Acts::ChargedPolicy' or "
                "'Acts::NeutralPolicy");

 public:
  /// Public typedefs
  /// Type of covariance matrix
  using CovMatrix_t = FreeSymMatrix;

  /// Construct track parameters for charged particles.
  ///
  /// @tparam T Type of the charge policy (ChargedPolicy)
  /// @param [in] cov The covariance matrix
  /// @param [in] parValues Vector with parameter values
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  SingleFreeTrackParameters(std::optional<CovMatrix_t> cov,
                            const FreeVector& parValues)
      : m_oParameters(std::move(cov), parValues),
        m_oChargePolicy(std::copysign(1., parValues[eFreeQOverP])) {}

  /// Construct track parameters for neutral particles.
  ///
  /// @tparam T Type of the charge policy (NeutralPolicy)
  /// @param [in] cov The covariance matrix
  /// @param [in] parValues Vector with parameter values
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  SingleFreeTrackParameters(std::optional<CovMatrix_t> cov,
                            const FreeVector& parValues)
      : m_oParameters(std::move(cov), parValues), m_oChargePolicy() {}

  // this class does not have a custom default constructor and thus should not
  // provide any custom default cstors, dstor, or assignment. see ISOCPP C.20.

  /// @brief Access all parameters
  ///
  /// @return Vector containing the store parameters
  FreeVector parameters() const { return m_oParameters.getParameters(); }

  /// @brief Access to a single parameter
  ///
  /// @tparam kIndex Identifier of the parameter index which will be retrieved
  ///
  /// @return Value of the requested parameter
  template <FreeParametersIndices kIndex,
            std::enable_if_t<kIndex<eFreeParametersSize, int> = 0> ParValue_t
                get() const {
    return m_oParameters.template getParameter<kIndex>();
  }

  /// @brief Access track parameter uncertainty
  ///
  /// @tparam par Identifier of the parameter uncertainty index which will
  /// be retrieved
  ///
  /// @return Value of the requested parameter uncertainty
  template <FreeParametersIndices kIndex,
            std::enable_if_t<kIndex<eFreeParametersSize, int> = 0> ParValue_t
                uncertainty() const {
    return m_oParameters.template getUncertainty<kIndex>();
  }

  /// @brief Access covariance matrix of track parameters
  ///
  /// @note The ownership of the covariance matrix is @b not transferred
  /// with this call.
  ///
  /// @return Raw pointer to covariance matrix (can be a nullptr)
  ///
  /// @sa ParameterSet::getCovariance
  const std::optional<CovMatrix_t>& covariance() const {
    return m_oParameters.getCovariance();
  }

  /// @brief access position in global coordinate system
  ///
  /// @return 3D vector with global position
  Vector3D position() const {
    return parameters().template segment<3>(eFreePos0);
  }

  /// @brief access momentum in global coordinate system
  ///
  /// @return 3D vector with global momentum
  Vector3D momentum() const {
    return parameters().template segment<3>(eFreeDir0) /
           std::abs(get<eFreeQOverP>());
  }

  /// @brief retrieve electric charge
  ///
  /// @return value of electric charge
  double charge() const { return m_oChargePolicy.getCharge(); }

  /// @brief retrieve time
  ///
  /// @return value of time
  double time() const { return get<eFreeTime>(); }

  /// @brief access to the internally stored FreeParameterSet
  ///
  /// @return FreeParameterSet object holding parameter values and their
  /// covariance matrix
  const FullFreeParameterSet& getParameterSet() const { return m_oParameters; }

  /// @brief Equality operator
  ///
  /// @param [in] rhs Object to compare `*this` to
  ///
  /// @return Boolean value whether the objects can be casted into each
  /// other and the content of the member variables is the same
  bool operator==(const SingleFreeTrackParameters& rhs) const {
    auto casted = dynamic_cast<decltype(this)>(&rhs);
    if (!casted) {
      return false;
    }

    return (m_oChargePolicy == casted->m_oChargePolicy &&
            m_oParameters == casted->m_oParameters);
  }

  /// @brief inequality operator
  ///
  /// @return `not (*this == rhs)`
  bool operator!=(const SingleFreeTrackParameters& rhs) const {
    return !(*this == rhs);
  }

  /// @brief Update of the parameterisation
  ///
  /// @tparam par The parameter index
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param newValue The new updaed value
  ///
  /// @note The context is not used here but makes the API consistent with
  /// @c SingleCurvilinearTrackParameters and @c SingleBoundTrackParameters
  template <FreeParametersIndices par,
            std::enable_if_t<par<eFreeParametersSize, int> = 0> void set(
                const GeometryContext& /*gctx*/, ParValue_t newValue) {
    m_oParameters.setParameter<par>(newValue);
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
                                  const SingleFreeTrackParameters& sfp) {
    sfp.print(out);
    return out;
  }

 private:
  FullFreeParameterSet
      m_oParameters;             ///< FreeParameterSet object holding the
                                 /// parameter values and covariance matrix
  ChargePolicy m_oChargePolicy;  ///< charge policy object distinguishing
                                 /// between charged and neutral tracks
};
}  // namespace Acts
// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/ParameterSet.hpp"
#include "Acts/EventData/detail/PrintParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include <iomanip>
#include <ostream>

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
  using Scalar = FreeParametersScalar;
  using ParametersVector = FreeVector;
  using CovarianceMatrix = FreeSymMatrix;

  /// Construct track parameters for charged particles.
  ///
  /// @tparam T Type of the charge policy (ChargedPolicy)
  /// @param [in] cov The covariance matrix
  /// @param [in] parValues Vector with parameter values
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  SingleFreeTrackParameters(std::optional<CovarianceMatrix> cov,
                            const ParametersVector& parValues)
      : m_oParameters(std::move(cov), parValues),
        m_oChargePolicy(std::copysign(1., parValues[eFreeQOverP])) {}

  /// Construct track parameters for neutral particles.
  ///
  /// @tparam T Type of the charge policy (NeutralPolicy)
  /// @param [in] cov The covariance matrix
  /// @param [in] parValues Vector with parameter values
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  SingleFreeTrackParameters(std::optional<CovarianceMatrix> cov,
                            const ParametersVector& parValues)
      : m_oParameters(std::move(cov), parValues), m_oChargePolicy() {}

  // this class does not have a custom default constructor and thus should not
  // provide any custom default cstors, dstor, or assignment. see ISOCPP C.20.

  /// @brief Access all parameters
  ///
  /// @return Vector containing the store parameters
  ParametersVector parameters() const { return m_oParameters.getParameters(); }

  /// @brief Access to a single parameter
  ///
  /// @tparam kIndex Identifier of the parameter index which will be retrieved
  ///
  /// @return Value of the requested parameter
  template <FreeParametersIndices kIndex>
  Scalar get() const {
    return m_oParameters.template getParameter<kIndex>();
  }

  /// @brief Access track parameter uncertainty
  ///
  /// @tparam kIndex Identifier of the uncertainty index which will be retrieved
  ///
  /// @return Value of the requested parameter uncertainty
  template <FreeParametersIndices kIndex>
  Scalar uncertainty() const {
    return m_oParameters.template getUncertainty<kIndex>();
  }

  /// @brief Access covariance matrix of track parameters
  ///
  /// @note The ownership of the covariance matrix is @b not transferred
  /// with this call.
  ///
  /// @sa ParameterSet::getCovariance
  const std::optional<CovarianceMatrix>& covariance() const {
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
  Scalar charge() const { return m_oChargePolicy.getCharge(); }

  /// @brief retrieve time
  ///
  /// @return value of time
  Scalar time() const { return get<eFreeTime>(); }

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
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param newValue The new updaed value
  /// @tparam kIndex Identifier of the parameter index which will be set
  ///
  /// @note The context is not used here but makes the API consistent with
  /// @c SingleCurvilinearTrackParameters and @c SingleBoundTrackParameters
  template <FreeParametersIndices kIndex>
  void set(const GeometryContext& /*gctx*/, Scalar newValue) {
    m_oParameters.setParameter<kIndex>(newValue);
  }

 private:
  FullFreeParameterSet
      m_oParameters;             ///< FreeParameterSet object holding the
                                 /// parameter values and covariance matrix
  ChargePolicy m_oChargePolicy;  ///< charge policy object distinguishing
                                 /// between charged and neutral tracks

  /// Print information to the output stream.
  friend std::ostream& operator<<(std::ostream& os,
                                  const SingleFreeTrackParameters& tp) {
    detail::printFreeParameters(
        os, tp.parameters(),
        tp.covariance().has_value() ? &tp.covariance().value() : nullptr);
    return os;
  }
};

}  // namespace Acts

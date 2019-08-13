// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <type_traits>
#include "Acts/EventData/FreeParametersBase.hpp"
#include "Acts/EventData/detail/coordinate_transformations.hpp"

namespace Acts {

/// @class SingleFreeParameters
///
/// @brief base class for a single set of track parameters
///
/// This class implements the interface for charged/neutral track parameters for
/// the case that it represents a single set of track parameters
/// (opposed to a list of different sets of track parameters as used by
/// e.g. GSF or multi-track fitters).
///
/// The track parameters and their uncertainty are defined in local reference
/// frame which depends on the associated surface of the track parameters.
///
/// @tparam ChargePolicy type for distinguishing charged and neutral
/// tracks/particles
///         (must be either ChargedPolicy or NeutralPolicy)
template <class ChargePolicy>
class SingleFreeParameters : public FreeParametersBase {
  static_assert(std::is_same<ChargePolicy, ChargedPolicy>::value or
                    std::is_same<ChargePolicy, NeutralPolicy>::value,
                "ChargePolicy must either be 'Acts::ChargedPolicy' or "
                "'Acts::NeutralPolicy");

 public:
  // public typedef's

  /// vector type for stored track parameters
  using ParVector_t = typename FreeParametersBase::ParVector_t;

  /// type of covariance matrix
  using CovMatrix_t = typename FreeParametersBase::CovMatrix_t;

  /// type for unique pointer to covariance matrix
  using CovPtr_t = std::unique_ptr<const CovMatrix_t>;

  /// @brief default virtual destructor
  ~SingleFreeParameters() override = default;

  /// @brief virtual constructor
  SingleFreeParameters<ChargePolicy>* clone() const override = 0;

  /// @copydoc FreeParametersBase::position
  ActsVectorD<3> position() const final { return m_vPosition; }

  /// @copydoc FreeParametersBase::momentum
  ActsVectorD<3> momentum() const final { return m_vMomentum; }

  /// @brief equality operator
  ///
  /// @return @c true of both objects have the same charge policy, parameter
  /// values, position and momentum, otherwise @c false
  bool operator==(const ParametersBase& rhs) const override {
    auto casted = dynamic_cast<decltype(this)>(&rhs);
    if (!casted) {
      return false;
    }

    return (m_oChargePolicy == casted->m_oChargePolicy &&
            m_oTime == casted->m_oTime &&
            m_oParameters == casted->m_oParameters &&
            m_vPosition == casted->m_vPosition &&
            m_vMomentum == casted->m_vMomentum);
  }

  /// @copydoc FreeParametersBase::charge
  double charge() const final { return m_oChargePolicy.getCharge(); }

  /// @copydoc FreeParametersBase::time
  double time() const final { return m_oTime; }

  /// @copydoc FreeParametersBase::getParameterSet
  const FullParameterSet& getParameterSet() const final {
    return m_oParameters;
  }

 protected:
  /// @brief standard constructor for track parameters of charged particles
  ///
  /// @param cov unique pointer to covariance matrix (nullptr is accepted)
  /// @param parValues vector with parameter values
  /// @param position 3D vector with global position
  /// @param momentum 3D vector with global momentum
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  SingleFreeParameters(CovPtr_t cov, const ParVector_t& parValues,
                        const ActsVectorD<3>& position,
                        const ActsVectorD<3>& momentum)
      : FreeParametersBase(),
        m_oChargePolicy(
            detail::coordinate_transformation::parameters2charge(parValues)),
        m_oTime(detail::coordinate_transformation::parameters2time(parValues)),
        m_oParameters(std::move(cov), parValues),
        m_vPosition(position),
        m_vMomentum(momentum) {}

  /// @brief standard constructor for track parameters of neutral particles
  ///
  /// @param cov unique pointer to covariance matrix (nullptr is accepted)
  /// @param parValues vector with parameter values
  /// @param position 3D vector with global position
  /// @param momentum 3D vector with global momentum
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  SingleFreeParameters(CovPtr_t cov, const ParVector_t& parValues,
                        const ActsVectorD<3>& position,
                        const ActsVectorD<3>& momentum)
      : FreeParametersBase(),
        m_oChargePolicy(),
        m_oTime(detail::coordinate_transformation::parameters2time(parValues)),
        m_oParameters(std::move(cov), parValues),
        m_vPosition(position),
        m_vMomentum(momentum) {}

  /// @brief default copy constructor
  SingleFreeParameters(const SingleFreeParameters<ChargePolicy>& copy) =
      default;

  /// @brief default move constructor
  SingleFreeParameters(SingleFreeParameters<ChargePolicy>&& copy) = default;

  /// @brief copy assignment operator
  ///
  /// @param rhs object to be copied
  SingleFreeParameters<ChargePolicy>& operator=(
      const SingleFreeParameters<ChargePolicy>& rhs) {
    // check for self-assignment
    if (this != &rhs) {
      m_oChargePolicy = rhs.m_oChargePolicy;
      m_oTime = rhs.m_oTime;
      m_oParameters = rhs.m_oParameters;
      m_vPosition = rhs.m_vPosition;
      m_vMomentum = rhs.m_vMomentum;
    }

    return *this;
  }

  /// @brief move assignment operator
  ///
  /// @param rhs object to be movied into `*this`
  SingleFreeParameters<ChargePolicy>& operator=(
      SingleFreeParameters<ChargePolicy>&& rhs) {
    // check for self-assignment
    if (this != &rhs) {
      m_oChargePolicy = std::move(rhs.m_oChargePolicy);
      m_oTime = std::move(rhs.m_oTime);
      m_oParameters = std::move(rhs.m_oParameters);
      m_vPosition = std::move(rhs.m_vPosition);
      m_vMomentum = std::move(rhs.m_vMomentum);
    }

    return *this;
  }

  /// @copydoc FreeParametersBase::getParameterSet
  virtual FullFreeParameterSet& getParameterSet() final { return m_oParameters; }

  /// @brief update global momentum from current parameter values
  ///
  ///
  /// @param[in] gctx is the Context object that is forwarded to the surface
  ///            for local to global coordinate transformation
  ///
  /// @note This function is triggered when called with an argument of a type
  ///       different from Acts::local_parameter
  template <typename T>
  void updateGlobalCoordinates(const GeometryContext& /*gctx*/,
                               const T& /*unused*/) {
    m_vMomentum = detail::coordinate_transformation::parameters2globalMomentum(
        getParameterSet().getParameters());
    m_oTime = detail::coordinate_transformation::parameters2time(
        getParameterSet().getParameters());
  }

  /// @brief update global position from current parameter values
  ///
  /// @note This function is triggered when called with an argument of a type
  /// Acts::local_parameter
  void updateGlobalCoordinates(const GeometryContext& gctx,
                               const local_parameter& /*unused*/) {
    m_vPosition = detail::coordinate_transformation::parameters2globalPosition(
        gctx, getParameterSet().getParameters(), this->referenceSurface());
  }

  ChargePolicy m_oChargePolicy;    ///< charge policy object distinguishing
                                   /// between charged and neutral tracks
  double m_oTime;                  ///< time of the track parametrisation
  FullFreeParameterSet m_oParameters;  ///< ParameterSet object holding the
                                   /// parameter values and covariance matrix
  FreeVector
  ActsVectorD<3> m_vPosition;      ///< 3D vector with global position
  ActsVectorD<3> m_vMomentum;      ///< 3D vector with global momentum
};
}  // namespace Acts
// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// STL include(s)
#include <ostream>

// Acts includes
#include "Acts/EventData/ParametersBase.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class TrackParametersBase
///
/// @brief base class for track parameters
///
/// This is a base class for neutral and charged track parameters.
/// The position and the momentum are both given in the global coordinate
/// system. The track parameters and their uncertainty are defined in local
/// reference frame which depends on the associated surface
/// of the track parameters.
template <class ChargePolicy>
class SingleFreeparameters : public ParametersBase {

static_assert(std::is_same<ChargePolicy, ChargedPolicy>::value or
                    std::is_same<ChargePolicy, NeutralPolicy>::value,
                "ChargePolicy must either be 'Acts::ChargedPolicy' or "
                "'Acts::NeutralPolicy");
 public:
  // public typedef's

  /// vector type for stored track parameters
  using ParVector_t = FreeVector;

  /// type of covariance matrix
  using CovMatrix_t = FreeSymMatrix;

  /// type for unique pointer to covariance matrix
  using CovPtr_t = std::unique_ptr<const CovMatrix_t>;
 
   /// @brief default virtual destructor
  ~SingleFreeparameters() override = default;

  /// @brief virtual constructor
  SingleFreeparameters<ChargePolicy>* clone() const override = 0;
   
  /// @brief access track parameters
  ///
  /// @return Eigen vector of dimension Acts::BoundParsDim with values of the
  /// track parameters
  ///         (in the order as defined by the ParID_t enumeration)
  ParVector_t parameters() const = 0;

  /// @brief access track parameter
  ///
  /// @tparam par identifier of track parameter which is to be retrieved
  ///
  /// @return value of the requested track parameter
  ///
  /// @sa ParameterSet::get
  template <unsigned int par>
  ParValue_t get() const  = 0;

  /// @brief access track parameter uncertainty
  ///
  /// @tparam par identifier of track parameter which is to be retrieved
  ///
  /// @return value of the requested track parameter uncertainty
  template <unsigned int par>
  ParValue_t uncertainty() const = 0;
  
  /// @brief access covariance matrix of track parameters
  ///
  /// @note The ownership of the covariance matrix is @b not transferred with
  /// this call.
  ///
  /// @return raw pointer to covariance matrix (can be a nullptr)
  ///
  /// @sa ParameterSet::getCovariance
  const CovMatrix_t* covariance() const = 0;




  /// @copydoc TrackParametersBase::position
  ActsVectorD<3> position() const final { return m_vPosition; }

  /// @copydoc TrackParametersBase::momentum
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

  /// @copydoc TrackParametersBase::charge
  double charge() const final { return m_oChargePolicy.getCharge(); }

  /// @copydoc TrackParametersBase::time
  double time() const final { return m_oTime; }




 protected:
  /// @brief standard constructor for track parameters of charged particles
  ///
  /// @param cov unique pointer to covariance matrix (nullptr is accepted)
  /// @param parValues vector with parameter values
  /// @param position 3D vector with global position
  /// @param momentum 3D vector with global momentum
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  SingleFreeparameters(CovPtr_t cov, const ParVector_t& parValues,
                        const ActsVectorD<3>& position,
                        const ActsVectorD<3>& momentum)
      : TrackParametersBase(),
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
  SingleFreeparameters(CovPtr_t cov, const ParVector_t& parValues,
                        const ActsVectorD<3>& position,
                        const ActsVectorD<3>& momentum)
      : TrackParametersBase(),
        m_oChargePolicy(),
        m_oTime(detail::coordinate_transformation::parameters2time(parValues)),
        m_oParameters(std::move(cov), parValues),
        m_vPosition(position),
        m_vMomentum(momentum) {}

  /// @brief default copy constructor
  SingleFreeparameters(const SingleTrackParameters<ChargePolicy>& copy) =
      default;

  /// @brief default move constructor
  SingleFreeparameters(SingleTrackParameters<ChargePolicy>&& copy) = default;

  /// @brief copy assignment operator
  ///
  /// @param rhs object to be copied
  SingleFreeparameters<ChargePolicy>& operator=(
      const SingleFreeparameters<ChargePolicy>& rhs) {
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
  SingleFreeparameters<ChargePolicy>& operator=(
      SingleFreeparameters<ChargePolicy>&& rhs) {
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



   
protected:
  /// @brief print information to output stream
  ///
  /// @return modified output stream object
  virtual std::ostream& print(std::ostream& sl) const;
  
  ChargePolicy m_oChargePolicy;    ///< charge policy object distinguishing
                                   /// between charged and neutral tracks                        
  ParVector_t m_parameters; ///< Parameter vector
  CovMatrix_t m_covariance; ///< Covariance matrix
  
};
}  // namespace Acts
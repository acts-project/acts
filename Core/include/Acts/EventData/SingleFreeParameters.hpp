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
#include <iomanip>

// Acts includes
#include "Acts/EventData/ParametersBase.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/EventData/detail/coordinate_transformations.hpp"

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
class SingleFreeParameters : public ParametersBase {

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
  ~SingleFreeParameters() override = default;

  /// @brief virtual constructor
  SingleFreeParameters<ChargePolicy>* clone() const override {
	return new SingleFreeParameters<ChargePolicy>(*this);
  } // TODO: This probably requires copying from curvilinear/bound, too
  
  /// @brief standard constructor for track parameters of charged particles
  ///
  /// @param cov unique pointer to covariance matrix (nullptr is accepted)
  /// @param parValues vector with parameter values
  /// @param position 3D vector with global position
  /// @param momentum 3D vector with global momentum
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  SingleFreeParameters(CovPtr_t cov, const ParVector_t& parValues)
      : ParametersBase(),
        m_oChargePolicy(
            detail::coordinate_transformation::parameters2charge(parValues)),
        m_covariance(std::move(cov)),
        m_parameters(parValues) {}

  /// @brief standard constructor for track parameters of neutral particles
  ///
  /// @param cov unique pointer to covariance matrix (nullptr is accepted)
  /// @param parValues vector with parameter values
  /// @param position 3D vector with global position
  /// @param momentum 3D vector with global momentum
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  SingleFreeParameters(CovPtr_t cov, const ParVector_t& parValues)
      : ParametersBase(),
        m_oChargePolicy(),
        m_covariance(std::move(cov)),
        m_parameters(parValues) {}

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
      m_parameters = rhs.m_parameters;
      m_covariance = (rhs.m_covariance
             ? std::make_unique<const CovMatrix_t>(*rhs.m_covariance)
             : nullptr);
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
      m_parameters = std::move(rhs.m_parameters);
      m_covariance = std::move(rhs.m_covariance);
    }

    return *this;
  }  
  
  /// @brief access track parameters
  ///
  /// @return Eigen vector of dimension Acts::BoundParsDim with values of the
  /// track parameters
  ///         (in the order as defined by the ParID_t enumeration)
  ParVector_t parameters() const { return m_parameters;}

  /// @brief access track parameter
  ///
  /// @tparam par identifier of track parameter which is to be retrieved
  ///
  /// @return value of the requested track parameter
  ///
  /// @sa ParameterSet::get
  template <unsigned int par,
            std::enable_if_t<par < FreeParsDim, int> = 0>
  ParValue_t get() const  { return m_parameters(par);}

  /// @brief access track parameter uncertainty
  ///
  /// @tparam par identifier of track parameter which is to be retrieved
  ///
  /// @return value of the requested track parameter uncertainty
  template <unsigned int par,
            std::enable_if_t<par < FreeParsDim, int> = 0>
  ParValue_t uncertainty() const { return m_covariance->coeff(par, par);}
  
  /// @brief access covariance matrix of track parameters
  ///
  /// @note The ownership of the covariance matrix is @b not transferred with
  /// this call.
  ///
  /// @return raw pointer to covariance matrix (can be a nullptr)
  ///
  /// @sa ParameterSet::getCovariance
  const CovMatrix_t* covariance() const { return m_covariance;}

  /// @copydoc TrackParametersBase::position
  ActsVectorD<3> position() const final { return m_parameters.template head<3>(); }

  /// @copydoc TrackParametersBase::momentum
  ActsVectorD<3> momentum() const final { return m_parameters.template segment<3>(4); } // TODO: this is the direction, not the momentum

  /// @brief equality operator
  ///
  /// @return @c true of both objects have the same charge policy, parameter
  /// values, position and momentum, otherwise @c false
  bool operator==(const ParametersBase& rhs) const override {
    auto casted = dynamic_cast<decltype(this)>(&rhs);
    if (!casted) {
      return false;
    }
	
	// both have covariance matrices set
    if ((m_covariance && casted->m_covariance) &&
        (*m_covariance != *casted->m_covariance)) {
      return false;
    }
    // only one has a covariance matrix set
    if ((m_covariance && !casted->m_covariance) ||
        (!m_covariance && casted->m_covariance)) {
      return false;
    }
    
    return (m_oChargePolicy == casted->m_oChargePolicy &&
            m_parameters == casted->m_parameters);
  }

  /// @copydoc TrackParametersBase::charge
  double charge() const final { return m_oChargePolicy.getCharge(); }

  /// @copydoc TrackParametersBase::time
  double time() const final { return m_parameters(3); }

  /// @brief update of the track parameterisation
  /// only possible on non-const objects, enable for local parameters
  ///
  /// @tparam ParID_t The parameter type
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param newValue The new updaed value
  ///
  /// For curvilinear parameters the local parameters are forced to be
  /// (0,0), hence an update is an effective shift of the reference
  template <unsigned int par,
            std::enable_if_t<par < FreeParsDim, int> = 0>
  void set(const GeometryContext& /*unused*/, ParValue_t newValue) {
	m_parameters(par) = newValue;
  }

  /// @brief print information to output stream
  ///
  /// @return modified output stream object
  std::ostream& print(std::ostream& sl) const {
	  // set stream output format
	  auto old_precision = sl.precision(7);
	  auto old_flags = sl.setf(std::ios::fixed);

	  sl << " * FreeTrackParameters: ";
	  sl << parameters().transpose() << std::endl;
	  sl << " * charge: " << charge() << std::endl;
	  if (covariance() != nullptr) {
		sl << " * covariance matrix:\n" << *covariance() << std::endl;
	  } else {
		sl << " * covariance matrix:\n" << covariance() << std::endl;
	  }

	  // reset stream format
	  sl.precision(old_precision);
	  sl.setf(old_flags);

	  return sl;
  }

private:
  ChargePolicy m_oChargePolicy;    ///< charge policy object distinguishing
                                   /// between charged and neutral tracks                        
  ParVector_t m_parameters; ///< Parameter vector // TODO: Could there be 3D references to position/direction?
  CovPtr_t m_covariance; ///< Covariance matrix
  Vector3D& m_position = m_parameters.template head<3>(); // TODO: test these two lines
  Vector3D& m_direction = m_parameters.template segment<3>(4);
};

}  // namespace Acts
// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <type_traits>

#include "Acts/EventData/ChargePolicy.hpp"
#include "Acts/EventData/ParameterSet.hpp"
#include "Acts/EventData/detail/coordinate_transformations.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class SingleTrackParameters
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
class SingleTrackParameters {
  static_assert(std::is_same<ChargePolicy, ChargedPolicy>::value or
                    std::is_same<ChargePolicy, NeutralPolicy>::value,
                "ChargePolicy must either be 'Acts::ChargedPolicy' or "
                "'Acts::NeutralPolicy");

 public:
  // public typedef's

  /// vector type for stored track parameters
  using ParVector_t = BoundVector;

  /// type of covariance matrix
  using CovMatrix_t = BoundSymMatrix;

  /// @brief default destructor
  virtual ~SingleTrackParameters() = default;

  /// @brief access position in global coordinate system
  ///
  /// @return 3D vector with global position
  Vector3D position() const { return m_vPosition; }

  /// @brief access momentum in global coordinate system
  ///
  /// @return 3D vector with global momentum
  Vector3D momentum() const { return m_vMomentum; }

  /// @brief equality operator
  ///
  /// @return @c true of both objects have the same charge policy, parameter
  /// values, position and momentum, otherwise @c false
  bool operator==(const SingleTrackParameters& rhs) const {
    auto casted = dynamic_cast<decltype(this)>(&rhs);
    if (!casted) {
      return false;
    }

    return (m_oChargePolicy == casted->m_oChargePolicy &&
            m_oParameters == casted->m_oParameters &&
            m_vPosition == casted->m_vPosition &&
            m_vMomentum == casted->m_vMomentum);
  }

  /// @brief retrieve electric charge
  ///
  /// @return value of electric charge
  double charge() const { return m_oChargePolicy.getCharge(); }

  /// @brief retrieve time
  ///
  /// @return value of time
  double time() const { return get<ParDef::eT>(); }

  /// @brief access to the internally stored ParameterSet
  ///
  /// @return ParameterSet object holding parameter values and their covariance
  /// matrix
  const FullParameterSet& getParameterSet() const { return m_oParameters; }

  /// @brief access associated surface defining the coordinate system for track
  ///        parameters and their covariance
  ///
  /// @return associated surface
  virtual const Surface& referenceSurface() const = 0;

  /// @brief access covariance matrix of track parameters
  ///
  /// @note The ownership of the covariance matrix is @b not transferred with
  /// this call.
  ///
  /// @return raw pointer to covariance matrix (can be a nullptr)
  ///
  /// @sa ParameterSet::getCovariance
  const std::optional<CovMatrix_t>& covariance() const {
    return getParameterSet().getCovariance();
  }

  /// @brief access track parameters
  ///
  /// @return Eigen vector of dimension Acts::eBoundParametersSize with values
  /// of the track parameters
  ///         (in the order as defined by the ParID_t enumeration)
  ParVector_t parameters() const { return getParameterSet().getParameters(); }

  /// @brief access track parameter
  ///
  /// @tparam par identifier of track parameter which is to be retrieved
  ///
  /// @return value of the requested track parameter
  ///
  /// @sa ParameterSet::get
  template <ParID_t par>
  ParValue_t get() const {
    return getParameterSet().template getParameter<par>();
  }

  /// @brief access track parameter uncertainty
  ///
  /// @tparam par identifier of track parameter which is to be retrieved
  ///
  /// @return value of the requested track parameter uncertainty
  template <ParID_t par>
  ParValue_t uncertainty() const {
    return getParameterSet().template getUncertainty<par>();
  }

  /// @brief convenience method to retrieve transverse momentum
  double pT() const { return VectorHelpers::perp(momentum()); }

  /// @brief convenience method to retrieve pseudorapidity
  double eta() const { return VectorHelpers::eta(momentum()); }

  FullParameterSet& getParameterSet() { return m_oParameters; }

  /// @brief output stream operator
  ///
  /// Prints information about this object to the output stream using the
  /// virtual
  /// TrackParameters::print method.
  ///
  /// @return modified output stream object
  friend std::ostream& operator<<(std::ostream& out,
                                  const SingleTrackParameters& stp) {
    stp.print(out);
    return out;
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
  SingleTrackParameters(std::optional<CovMatrix_t> cov,
                        const ParVector_t& parValues, const Vector3D& position,
                        const Vector3D& momentum)
      : m_oChargePolicy(
            detail::coordinate_transformation::parameters2charge(parValues)),
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
  SingleTrackParameters(std::optional<CovMatrix_t> cov,
                        const ParVector_t& parValues, const Vector3D& position,
                        const Vector3D& momentum)
      : m_oChargePolicy(),
        m_oParameters(std::move(cov), parValues),
        m_vPosition(position),
        m_vMomentum(momentum) {}

  /// @brief default copy constructor
  SingleTrackParameters(const SingleTrackParameters<ChargePolicy>& copy) =
      default;

  /// @brief default move constructor
  SingleTrackParameters(SingleTrackParameters<ChargePolicy>&& copy) = default;

  /// @brief copy assignment operator
  ///
  /// @param rhs object to be copied
  SingleTrackParameters<ChargePolicy>& operator=(
      const SingleTrackParameters<ChargePolicy>& rhs) {
    // check for self-assignment
    if (this != &rhs) {
      m_oChargePolicy = rhs.m_oChargePolicy;
      m_oParameters = rhs.m_oParameters;
      m_vPosition = rhs.m_vPosition;
      m_vMomentum = rhs.m_vMomentum;
    }

    return *this;
  }

  /// @brief move assignment operator
  ///
  /// @param rhs object to be movied into `*this`
  SingleTrackParameters<ChargePolicy>& operator=(
      SingleTrackParameters<ChargePolicy>&& rhs) {
    // check for self-assignment
    if (this != &rhs) {
      m_oChargePolicy = std::move(rhs.m_oChargePolicy);
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

  /// @brief print information to output stream
  ///
  /// @return modified output stream object
  std::ostream& print(std::ostream& sl) const {
    // set stream output format
    auto old_precision = sl.precision(7);
    auto old_flags = sl.setf(std::ios::fixed);

    sl << " * TrackParameters: ";
    sl << parameters().transpose() << std::endl;
    sl << " * charge: " << charge() << std::endl;
    if (covariance()) {
      sl << " * covariance matrix:\n" << *covariance() << std::endl;
    } else {
      sl << " * covariance matrix:\nnull" << std::endl;
    }
    sl << " * corresponding global parameters:" << std::endl;
    sl << " *    position  (x y z) = (" << position().transpose() << ")"
       << std::endl;
    sl << " *    momentum  (px py pz) = (" << momentum().transpose() << ")"
       << std::endl;

    // reset stream format
    sl.precision(old_precision);
    sl.setf(old_flags);

    return sl;
  }

  ChargePolicy m_oChargePolicy;    ///< charge policy object distinguishing
                                   /// between charged and neutral tracks
  FullParameterSet m_oParameters;  ///< ParameterSet object holding the
                                   /// parameter values and covariance matrix
  Vector3D m_vPosition;            ///< 3D vector with global position
  Vector3D m_vMomentum;            ///< 3D vector with global momentum
};

}  // namespace Acts

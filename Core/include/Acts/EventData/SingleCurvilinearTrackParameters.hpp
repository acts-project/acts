// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <memory>
#include "Acts/EventData/SingleTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

namespace Acts {

/// @brief Charged and Neutrial Curvilinear Track representation
/// This is a single-component representation
///
/// @note the curvilinear representation is bound to the curvilinear
/// planar surface represenation. I.e. the local parameters are by
/// construction (0,0), the curvilinear surface is characterised by
/// being perpendicular to the track direction. It's internal frame
/// is constructed with the help of the global z axis.
template <typename ChargePolicy>
class SingleCurvilinearTrackParameters
    : public SingleTrackParameters<ChargePolicy> {
 public:
  /// type of covariance matrix
  using CovPtr_t = typename SingleTrackParameters<ChargePolicy>::CovPtr_t;

  /// @brief constructor for curvilienear representation
  /// This is the constructor from global parameters, enabled only
  /// for charged representations.
  ///
  /// @param[in] cov The covariance matrix w.r.t. curvilinear frame
  /// @param[in] position The global position of this track parameterisation
  /// @param[in] momentum The global momentum of this track parameterisation
  /// @param[in] dCharge The charge of this track parameterisation
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  SingleCurvilinearTrackParameters(CovPtr_t cov, const ActsVectorD<3>& position,
                                   const ActsVectorD<3>& momentum,
                                   double dCharge, double dTime)
      : SingleTrackParameters<ChargePolicy>(
            std::move(cov),
            detail::coordinate_transformation::global2curvilinear(
                position, momentum, dCharge, dTime),
            position, momentum),
        m_upSurface(Surface::makeShared<PlaneSurface>(position, momentum)) {}

  /// @brief constructor for curvilienear representation
  /// This is the constructor from global parameters, enabled only
  /// for charged representations.
  ///
  /// @param[in] cov The covariance matrix w.r.t. curvilinear frame
  /// @param[in] position The global position of this track parameterisation
  /// @param[in] momentum The global momentum of this track parameterisation
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  SingleCurvilinearTrackParameters(CovPtr_t cov, const ActsVectorD<3>& position,
                                   const ActsVectorD<3>& momentum, double dTime)
      : SingleTrackParameters<ChargePolicy>(
            std::move(cov),
            detail::coordinate_transformation::global2curvilinear(
                position, momentum, 0, dTime),
            position, momentum),
        m_upSurface(Surface::makeShared<PlaneSurface>(position, momentum)) {}

  /// @brief copy constructor - charged/neutral
  /// @param[in] copy The source parameters
  SingleCurvilinearTrackParameters(
      const SingleCurvilinearTrackParameters<ChargePolicy>& copy)
      : SingleTrackParameters<ChargePolicy>(copy),
        m_upSurface(copy.m_upSurface)  // copy shared ptr
  {}

  /// @brief move constructor - charged/neutral
  /// @param[in] other The source parameters
  SingleCurvilinearTrackParameters(
      SingleCurvilinearTrackParameters<ChargePolicy>&& other)
      : SingleTrackParameters<ChargePolicy>(std::move(other)),
        m_upSurface(std::move(other.m_upSurface)) {}

  /// @brief desctructor
  ~SingleCurvilinearTrackParameters() override = default;

  /// @brief copy assignment operator - charged/netural
  /// virtual constructor for type creation without casting
  SingleCurvilinearTrackParameters<ChargePolicy>& operator=(
      const SingleCurvilinearTrackParameters<ChargePolicy>& rhs) {
    // check for self-assignment
    if (this != &rhs) {
      SingleTrackParameters<ChargePolicy>::operator=(rhs);
      m_upSurface =
          Surface::makeShared<PlaneSurface>(this->position(), this->momentum());
    }
    return *this;
  }

  /// @brief move assignment operator - charged/netural
  /// virtual constructor for type creation without casting
  SingleCurvilinearTrackParameters<ChargePolicy>& operator=(
      SingleCurvilinearTrackParameters<ChargePolicy>&& rhs) {
    // check for self-assignment
    if (this != &rhs) {
      SingleTrackParameters<ChargePolicy>::operator=(std::move(rhs));
      m_upSurface = std::move(rhs.m_upSurface);
    }
    return *this;
  }

  /// @brief clone - charged/netural
  /// virtual constructor for type creation without casting
  SingleTrackParameters<ChargePolicy>* clone() const override {
    return new SingleCurvilinearTrackParameters<ChargePolicy>(*this);
  }

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
  template <ParID_t par,
            std::enable_if_t<std::is_same<typename par_type<par>::type,
                                          local_parameter>::value,
                             int> = 0>
  void set(const GeometryContext& gctx, ParValue_t newValue) {
    // set the parameter & update the new global position
    this->getParameterSet().template setParameter<par>(newValue);
    this->updateGlobalCoordinates(gctx, typename par_type<par>::type());
    // recreate the surface
    m_upSurface = Surface::makeShared<PlaneSurface>(
        this->position(), this->momentum().normalized());
    // reset to (0,0)
    this->getParameterSet().template setParameter<par>(0.);
  }

  /// @brief update of the track parameterisation
  /// only possible on non-const objects
  /// enable for parameters that are not local parameters
  /// @tparam ParID_t The parameter type
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param newValue The new updaed value
  ///
  /// For curvilinear parameters the directional change of parameters
  /// causes a recalculation of the surface
  template <ParID_t par,
            std::enable_if_t<not std::is_same<typename par_type<par>::type,
                                              local_parameter>::value,
                             int> = 0>
  void set(const GeometryContext& gctx, ParValue_t newValue) {
    this->getParameterSet().template setParameter<par>(newValue);
    this->updateGlobalCoordinates(gctx, typename par_type<par>::type());
    // recreate the surface
    m_upSurface = Surface::makeShared<PlaneSurface>(
        this->position(), this->momentum().normalized());
  }

  /// @brief access to the reference surface
  const Surface& referenceSurface() const final { return *m_upSurface; }

  /// @brief access to the measurement frame, i.e. the rotation matrix with
  /// respect to the global coordinate system, in which the local error
  /// is described.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///             It is ignored for Curvilinear parameters
  ///
  /// @note For a curvilinear track parameterisation this is identical to
  /// the rotation matrix of the intrinsic planar surface.
  RotationMatrix3D referenceFrame(const GeometryContext& gctx) const final {
    return m_upSurface->transform(gctx).linear();
  }

 private:
  std::shared_ptr<PlaneSurface> m_upSurface;
};
}  // namespace Acts
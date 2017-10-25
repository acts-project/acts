// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_SINGLECURVILINEARTRACKPARAMETERS_H
#define ACTS_SINGLECURVILINEARTRACKPARAMETERS_H

#include <memory>
#include "ACTS/EventData/SingleTrackParameters.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"

namespace Acts {

/// @brief Charged and Neutrial Curvilinear Track representation
/// This is a single-component representation
///
/// @note the curvilinear representation is bound to the curvilinear
/// planar surface represenation. I.e. the local parameters are by
/// construction (0,0), the curvilinear surface is characterised by
/// being perpenticular to the track direction. It's internal frame
/// is constructed with the help of the global z axis.
template <typename ChargePolicy>
class SingleCurvilinearTrackParameters
    : public SingleTrackParameters<ChargePolicy>
{
public:
  typedef typename SingleTrackParameters<ChargePolicy>::CovPtr_t
      CovPtr_t;  ///< type of covariance matrix

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
  SingleCurvilinearTrackParameters(CovPtr_t              cov,
                                   const ActsVectorD<3>& position,
                                   const ActsVectorD<3>& momentum,
                                   double                dCharge)
    : SingleTrackParameters<ChargePolicy>(
          std::move(cov),
          detail::coordinate_transformation::global2curvilinear(position,
                                                                momentum,
                                                                dCharge),
          position,
          momentum)
    , m_upSurface(new PlaneSurface(position, momentum))
  {
  }

  /// @brief constructor for curvilienear representation
  /// This is the constructor from global parameters, enabled only
  /// for charged representations.
  ///
  /// @param[in] cov The covariance matrix w.r.t. curvilinear frame
  /// @param[in] position The global position of this track parameterisation
  /// @param[in] momentum The global momentum of this track parameterisation
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  SingleCurvilinearTrackParameters(CovPtr_t              cov,
                                   const ActsVectorD<3>& position,
                                   const ActsVectorD<3>& momentum)
    : SingleTrackParameters<ChargePolicy>(
          std::move(cov),
          detail::coordinate_transformation::global2curvilinear(position,
                                                                momentum,
                                                                0),
          position,
          momentum)
    , m_upSurface(new PlaneSurface(position, momentum))
  {
  }

  /// @brief copy constructor - charged/neutral
  /// @param[in] copy The source parameters
  SingleCurvilinearTrackParameters(
      const SingleCurvilinearTrackParameters<ChargePolicy>& copy)
    : SingleTrackParameters<ChargePolicy>(copy)
    , m_upSurface(new PlaneSurface(this->position(), this->momentum()))
  {
  }

  /// @brief move constructor - charged/neutral
  /// @param[in] copy The source parameters
  SingleCurvilinearTrackParameters(
      SingleCurvilinearTrackParameters<ChargePolicy>&& copy)
    : SingleTrackParameters<ChargePolicy>(std::move(copy))
    , m_upSurface(std::move(copy.m_upSurface))
  {
  }

  /// @brief desctructor
  virtual ~SingleCurvilinearTrackParameters() = default;

  /// @brief copy assignment operator - charged/netural
  /// virtual constructor for type creation without casting
  SingleCurvilinearTrackParameters<ChargePolicy>&
  operator=(const SingleCurvilinearTrackParameters<ChargePolicy>& rhs)
  {
    // check for self-assignment
    if (this != &rhs) {
      SingleTrackParameters<ChargePolicy>::operator=(rhs);
      m_upSurface.reset(new PlaneSurface(this->position(), this->momentum()));
    }

    return *this;
  }

  /// @brief move assignment operator - charged/netural
  /// virtual constructor for type creation without casting
  SingleCurvilinearTrackParameters<ChargePolicy>&
  operator=(SingleCurvilinearTrackParameters<ChargePolicy>&& rhs)
  {
    // check for self-assignment
    if (this != &rhs) {
      SingleTrackParameters<ChargePolicy>::operator=(std::move(rhs));
      m_upSurface                                  = std::move(rhs.m_upSurface);
    }

    return *this;
  }

  /// @brief clone - charged/netural
  /// virtual constructor for type creation without casting
  virtual SingleTrackParameters<ChargePolicy>*
  clone() const override
  {
    return new SingleCurvilinearTrackParameters<ChargePolicy>(*this);
  }

  template <ParID_t par,
            std::enable_if_t<not std::is_same<typename par_type<par>::type,
                                              local_parameter>::value,
                             int> = 0>
  void
  set(ParValue_t newValue)
  {
    this->getParameterSet().template setParameter<par>(newValue);
    this->updateGlobalCoordinates(typename par_type<par>::type());
  }

  virtual const Surface&
  referenceSurface() const final override
  {
    return *m_upSurface;
  }

private:
  std::unique_ptr<const PlaneSurface> m_upSurface;
};
}  // end of namespace Acts

#endif  // ACTS_SINGLECURVILINEARTRACKPARAMETERS_H

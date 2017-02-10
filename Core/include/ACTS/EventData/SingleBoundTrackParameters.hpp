// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_SINGLEBOUNDTRACKPARAMETERS_H
#define ACTS_SINGLEBOUNDTRACKPARAMETERS_H 1

// ACTS includes
#include "ACTS/EventData/SingleTrackParameters.hpp"
#include "ACTS/Surfaces/Surface.hpp"

namespace Acts {
template <class ChargePolicy>
class SingleBoundTrackParameters : public SingleTrackParameters<ChargePolicy>
{
public:
  typedef typename SingleTrackParameters<ChargePolicy>::ParVector_t ParVector_t;
  typedef typename SingleTrackParameters<ChargePolicy>::CovPtr_t    CovPtr_t;

  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  SingleBoundTrackParameters(CovPtr_t           cov,
                             const ParVector_t& parValues,
                             const Surface&     surface)
    : SingleTrackParameters<ChargePolicy>(
          std::move(cov),
          parValues,
          detail::coordinate_transformation::parameters2globalPosition(
              parValues,
              surface),
          detail::coordinate_transformation::parameters2globalMomentum(
              parValues))
    , m_pSurface(surface.isFree() ? surface.clone() : &surface)
  {
  }

  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  SingleBoundTrackParameters(CovPtr_t              cov,
                             const ActsVectorD<3>& position,
                             const ActsVectorD<3>& momentum,
                             double                dCharge,
                             const Surface&        surface)
    : SingleTrackParameters<ChargePolicy>(
          std::move(cov),
          detail::coordinate_transformation::global2parameters(position,
                                                               momentum,
                                                               dCharge,
                                                               surface),
          position,
          momentum)
    , m_pSurface(surface.isFree() ? surface.clone() : &surface)
  {
  }

  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  SingleBoundTrackParameters(CovPtr_t           cov,
                             const ParVector_t& parValues,
                             const Surface&     surface)
    : SingleTrackParameters<ChargePolicy>(
          std::move(cov),
          parValues,
          detail::coordinate_transformation::parameters2globalPosition(
              parValues,
              surface),
          detail::coordinate_transformation::parameters2globalMomentum(
              parValues))
    , m_pSurface(surface.isFree() ? surface.clone() : &surface)
  {
  }

  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  SingleBoundTrackParameters(CovPtr_t              cov,
                             const ActsVectorD<3>& position,
                             const ActsVectorD<3>& momentum,
                             const Surface&        surface)
    : SingleTrackParameters<ChargePolicy>(
          std::move(cov),
          detail::coordinate_transformation::global2parameters(position,
                                                               momentum,
                                                               0,
                                                               surface),
          position,
          momentum)
    , m_pSurface(surface.isFree() ? surface.clone() : &surface)
  {
  }

  /**
   * @brief copy constructor
   */
  SingleBoundTrackParameters(
      const SingleBoundTrackParameters<ChargePolicy>& copy)
    : SingleTrackParameters<ChargePolicy>(copy)
    , m_pSurface(copy.m_pSurface->isFree() ? copy.m_pSurface->clone()
                                           : copy.m_pSurface)
  {
  }

  /**
   * @brief move constructor
   */
  SingleBoundTrackParameters(SingleBoundTrackParameters<ChargePolicy>&& copy)
    : SingleTrackParameters<ChargePolicy>(std::move(copy))
    , m_pSurface(copy.m_pSurface)
  {
    copy.m_pSurface = 0;
  }

  virtual ~SingleBoundTrackParameters()
  {
    if (m_pSurface && m_pSurface->isFree()) delete m_pSurface;
  }

  /**
   * @brief copy assignment operator
   */
  SingleBoundTrackParameters<ChargePolicy>&
  operator=(const SingleBoundTrackParameters<ChargePolicy>& rhs)
  {
    // check for self-assignment
    if (this != &rhs) {
      SingleTrackParameters<ChargePolicy>::operator=(rhs);

      if (m_pSurface->isFree()) delete m_pSurface;

      m_pSurface
          = rhs.m_pSurface->isFree() ? rhs.m_pSurface->clone() : rhs.m_pSurface;
    }

    return *this;
  }

  /**
   * @brief move assignment operator
   */
  SingleBoundTrackParameters<ChargePolicy>&
  operator=(SingleBoundTrackParameters<ChargePolicy>&& rhs)
  {
    // check for self-assignment
    if (this != &rhs) {
      SingleTrackParameters<ChargePolicy>::operator=(std::move(rhs));

      if (m_pSurface->isFree()) delete m_pSurface;

      m_pSurface     = rhs.m_pSurface;
      rhs.m_pSurface = 0;
    }

    return *this;
  }

  virtual SingleBoundTrackParameters<ChargePolicy>*
  clone() const override
  {
    return new SingleBoundTrackParameters<ChargePolicy>(*this);
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
  referenceSurface() const final
  {
    return *m_pSurface;
  }

private:
  const Surface* m_pSurface;
};

}  // end of namespace Acts

#endif  // ACTS_SINGLEBOUNDTRACKPARAMETERS_H

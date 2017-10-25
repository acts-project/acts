// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_SINGLEBOUNDTRACKPARAMETERS_H
#define ACTS_SINGLEBOUNDTRACKPARAMETERS_H

#include "ACTS/EventData/SingleTrackParameters.hpp"
#include "ACTS/Surfaces/Surface.hpp"

namespace Acts {

/// @brief Charged and Neutrial Track Parameterisation classes bound to
/// to a reference surface. This is a single-component representation
///
/// @note surfaces are dealt as plain pointers as they might represent
/// actual detector elements with their respective changeing alignment,
/// this might be reviewed at later stage in the context of data locality.
/// Currently, the following strategy is implemented:
/// * surfaces with a lifetime longer than an event, i.e. detector elements,
///   or representing surfaces of the detector geometry are copied as plain
///   pointers and hence not deleted as their ownership staus with the geometry
/// * (free) surfaces that have no registered ownership are cloned
template <class ChargePolicy>
class SingleBoundTrackParameters : public SingleTrackParameters<ChargePolicy>
{
public:
  typedef typename SingleTrackParameters<ChargePolicy>::ParVector_t ParVector_t;
  typedef typename SingleTrackParameters<ChargePolicy>::CovPtr_t    CovPtr_t;

  /// @brief Constructor of track parameters bound to a surface
  /// This is the constructor from global parameters, enabled only
  /// for charged representations.
  ///
  /// The transformations declared int he coordinate_transformation
  /// yield the global parameters and momentum representation
  /// @param[in] cov The covaraniance matrix (optional, can be nullptr)
  ///            it is given in the measurement frame
  /// @param[in] parValues The parameter vector
  /// @param[in] surface The reference surface the parameters are bound to
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

  /// @brief Constructor of track parameters bound to a surface
  /// This is the constructor from global parameters, enabled only
  /// for charged representations.
  ///
  /// The transformations declared int he coordinate_transformation
  /// yield the local parameters
  ///
  /// @param[in] cov The covaraniance matrix (optional, can be nullptr)
  ///            it is given in the curvilinear frame
  /// @param[in] position The global position of the track parameterisation
  /// @param[in] momentum The global momentum of the track parameterisation
  /// @param[in] dCharge The charge of the particle track parameterisation
  /// @param[in] surface The reference surface the parameters are bound to
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

  /// @brief Constructor of track parameters bound to a surface
  /// This is the constructor from global parameters, enabled only
  /// for neutral representations.
  ///
  /// The transformations declared int he coordinate_transformation
  /// yield the global parameters and momentum representation
  /// @param[in] cov The covaraniance matrix (optional, can be nullptr)
  ///            it is given in the measurement frame
  /// @param[in] parValues The parameter vector
  /// @param[in] surface The reference surface the parameters are bound to
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

  /// @brief Constructor of track parameters bound to a surface
  /// This is the constructor from global parameters, enabled only
  /// for neutral representations.
  ///
  /// The transformations declared int he coordinate_transformation
  /// yield the local parameters
  ///
  /// @param[in] cov The covaraniance matrix (optional, can be nullptr)
  ///            it is given in the curvilinear frame
  /// @param[in] position The global position of the track parameterisation
  /// @param[in] momentum The global momentum of the track parameterisation
  /// @param[in] dCharge The charge of the particle track parameterisation
  /// @param[in] surface The reference surface the parameters are bound to
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

  /// @brief copy constructor  - charged/neutral
  /// @param[in] copy The source parameters
  SingleBoundTrackParameters(
      const SingleBoundTrackParameters<ChargePolicy>& copy)
    : SingleTrackParameters<ChargePolicy>(copy)
    , m_pSurface(copy.m_pSurface->isFree() ? copy.m_pSurface->clone()
                                           : copy.m_pSurface)
  {
  }

  /// @brief move constructor - charged/neutral
  /// @param[in] copy The source parameters
  SingleBoundTrackParameters(SingleBoundTrackParameters<ChargePolicy>&& copy)
    : SingleTrackParameters<ChargePolicy>(std::move(copy))
    , m_pSurface(copy.m_pSurface)
  {
    copy.m_pSurface = 0;
  }

  /// @brief desctructor - charged/neutral
  /// checks if the surface is free and in such a case deletes it
  virtual ~SingleBoundTrackParameters()
  {
    if (m_pSurface && m_pSurface->isFree()) delete m_pSurface;
  }

  /// @brief copy assignment operator - charged/neutral
  /// checks if the surface is free and in such a case delete-clones it
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

  /// @brief move assignment operator - charged/neutral
  /// checks if the surface is free and in such a case delete-clones it
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

  /// @brief clone - charged/netural
  /// virtual constructor for type creation without casting
  virtual SingleBoundTrackParameters<ChargePolicy>*
  clone() const override
  {
    return new SingleBoundTrackParameters<ChargePolicy>(*this);
  }

  /// @brief set method for parameter updates
  /// enabled only if the parameter set is the same
  /// obviously only allowed on non-cost objects
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

  /// @brief access method to the reference surface
  virtual const Surface&
  referenceSurface() const final override
  {
    return *m_pSurface;
  }

  /// @brief access to the reference frame
  /// this is the logical frame where the covariance matrix
  /// is bound to - for curvilinear (== planar) surfaces
  /// this is given by the rotation matrix of the plane transform
  virtual RotationMatrix3D
  referenceFrame() const final override
  {
    /// @todo fix & update
    return RotationMatrix3D::Zero();
  }

private:
  const Surface* m_pSurface;
};

}  // end of namespace Acts

#endif  // ACTS_SINGLEBOUNDTRACKPARAMETERS_H

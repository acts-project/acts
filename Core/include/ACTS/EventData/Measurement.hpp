// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_MEASUREMENT_H
#define ACTS_MEASUREMENT_H 1

// STL include(s)
#include <memory>
#include <ostream>
#include <type_traits>
#include <utility>

// ACTS includes
#include "ACTS/EventData/ParameterSet.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/EventData/detail/fittable_type_generator.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"

namespace Acts {
// forward declarations
class Surface;

/**
 * @brief base class for Measurements
 *
 * This class describes the measurement of track parameters at a certain Surface
 * in the
 * TrackingGeometry.
 *
 * @note Identifier must be copy-constructible, move-constructible,
 * copy-assignable and move-assignable.
 *
 * @test The behavior of this class is tested in the following unit test:
 *       - \link Acts::Test::BOOST_AUTO_TEST_CASE(measurement_initialization)
 * initialization\endlink
 *
 * @tparam Identifier identification object for this measurement
 * @tparam params     parameter pack containing the measured parameters
 */
template <typename Identifier, ParID_t... params>
class Measurement
{
  // check type conditions
  static_assert(std::is_copy_constructible<Identifier>::value,
                "'Identifier' must be copy-constructible");
  static_assert(std::is_move_constructible<Identifier>::value,
                "'Identifier' must be move-constructible");
  static_assert(std::is_copy_assignable<Identifier>::value,
                "'Identifier' must be copy-assignable");
  static_assert(std::is_move_assignable<Identifier>::value,
                "'Identifier' must be move-assignable");

private:
  // private typedef's
  typedef ParameterSet<params...>
      ParSet_t;  ///< type of the underlying ParameterSet object

public:
  /// type of the vector containing the parameter values
  typedef typename ParSet_t::ParVector_t ParVector_t;
  /// type of the covariance matrix of the measurement
  typedef typename ParSet_t::CovMatrix_t CovMatrix_t;
  /// matrix type for projecting full parameter vector onto local parameter
  /// space
  typedef typename ParSet_t::Projection_t Projection_t;

  /**
   * @brief standard constructor
   *
   * Interface class for all possible measurements.
   *
   * @note Only a reference to the given surface is stored. The user must
   * ensure
   * that the lifetime of the @c Surface
   *       object surpasses the lifetime of this Measurement object.<br />
   *       The given parameter values are interpreted as values to the
   * parameters as defined in the class template
   *       argument @c params.
   *
   * @attention The current design will fail if the in-memory location of
   * the @c
   * Surface object is changed (e.g.
   *            if it is stored in a container and this gets relocated).
   *
   * @param surface surface at which the measurement took place
   * @param id identification object for this measurement
   * @param cov covariance matrix of the measurement.
   * @param head,values consistent number of parameter values of the
   * measurement
   */
  template <typename... Tail>
  Measurement(const Surface&    surface,
              const Identifier& id,
              CovMatrix_t       cov,
              typename std::enable_if<sizeof...(Tail) + 1 == sizeof...(params),
                                      ParValue_t>::type head,
              Tail... values)
    : m_oParameters(std::make_unique<CovMatrix_t>(std::move(cov)),
                    head,
                    values...)
    , m_pSurface(surface.isFree() ? const_cast<const Surface*>(surface.clone())
                                  : &surface)
    , m_oIdentifier(id)
  {
  }

  /**
   * @brief virtual destructor
   */
  virtual ~Measurement()
  {
    if (m_pSurface && m_pSurface->isFree()) {
      delete m_pSurface;
      m_pSurface = 0;
    }
  }

  /**
   * @brief copy constructor
   */
  Measurement(const Measurement<Identifier, params...>& copy)
    : m_oParameters(copy.m_oParameters)
    , m_pSurface(copy.m_pSurface->isFree()
                     ? const_cast<const Surface*>(copy.m_pSurface->clone())
                     : copy.m_pSurface)
    , m_oIdentifier(copy.m_oIdentifier)
  {
  }

  /**
   * @brief move constructor
   */
  Measurement(Measurement<Identifier, params...>&& rhs)
    : m_oParameters(std::move(rhs.m_oParameters))
    , m_pSurface(rhs.m_pSurface)
    , m_oIdentifier(std::move(rhs.m_oIdentifier))
  {
    rhs.m_pSurface = 0;
  }

  /**
   * @brief copy assignment operator
   */
  Measurement<Identifier, params...>&
  operator=(const Measurement<Identifier, params...>& rhs)
  {
    m_oParameters = rhs.m_oParameters;
    m_pSurface    = rhs.m_pSurface;
    m_oIdentifier = rhs.m_oIdentifier;

    return *this;
  }

  /**
   * @brief move assignment operator
   */
  Measurement<Identifier, params...>&
  operator=(Measurement<Identifier, params...>&& rhs)
  {
    m_oParameters = std::move(rhs.m_oParameters);
    m_pSurface    = rhs.m_pSurface;
    m_oIdentifier = std::move(rhs.m_oIdentifier);

    return *this;
  }

  /**
   * @brief retrieve stored value for given parameter
   *
   * @tparam parameter identifier for the parameter to be retrieved
   * @remark @c parameter must be part of the template parameter pack @c params.
   * Otherwise a compile-time
   *         error is generated.
   *
   * @return value of the stored parameter
   */
  template <ParID_t parameter>
  ParValue_t
  get() const
  {
    return m_oParameters.template getParameter<parameter>();
  }

  /**
   * @brief access vector with measured parameter values
   *
   * @return column vector whose size is equal to the dimensionality of this
   * Measurement. The values are
   *         given for the measured parameters in the order defined by the class
   * template argument @c params.
   */
  ParVector_t
  parameters() const
  {
    return m_oParameters.getParameters();
  }

  /**
   * @brief access covariance matrix of the measured parameter values
   *
   * @return covariance matrix of the measurement
   */
  CovMatrix_t
  covariance() const
  {
    return *m_oParameters.getCovariance();
  }

  /**
   * @brief retrieve stored uncertainty for given parameter
   *
   * @tparam parameter identifier for the parameter to be retrieved
   * @remark @c parameter must be part of the template parameter pack @c params.
   * Otherwise a compile-time
   *         error is generated.
   *
   * @return uncertainty \f$\sigma \ge 0\f$ for given parameter
   */
  template <ParID_t parameter>
  ParValue_t
  uncertainty() const
  {
    return m_oParameters.template getUncertainty<parameter>();
  }

  /**
   * @brief number of measured parameters
   *
   * @return number of measured parameters
   */
  static constexpr unsigned int
  size()
  {
    return ParSet_t::size();
  }

  /**
   * @brief access associated surface
   *
   * @pre The @c Surface object used to construct this @c Measurement object
   * must still be valid
   *      at the same memory location.
   *
   * @return reference to surface at which the measurement took place
   */
  const Acts::Surface&
  associatedSurface() const
  {
    return *m_pSurface;
  }

  /**
   * @brief access to global measurement identifier
   *
   * @return identifier object
   */
  Identifier
  identifier() const
  {
    return m_oIdentifier;
  }

  /**
   * @brief calculate residual with respect to given track parameters
   *
   * @note It is checked that the residual for non-local parameters are in valid
   * range (e.g.
   *       residuals in \f$\phi\f$ are corrected).
   *
   * @todo Implement check that TrackParameters are defined at the same Surface
   * as the Measurement is.
   * @todo Implement validity check for residuals of local parameters.
   *
   * @param trackPars reference TrackParameters object
   *
   * @return vector with the residual parameter values (in valid range)
   *
   * @sa ParameterSet::residual
   */
  ParVector_t
  residual(const TrackParameters& trackPars) const
  {
    return m_oParameters.residual(trackPars.getParameterSet());
  }

  /**
   * @brief equality operator
   *
   * @return @c true if parameter sets and associated surfaces compare equal,
   * otherwise @c false
   */
  virtual bool
  operator==(const Measurement<Identifier, params...>& rhs) const
  {
    return ((m_oParameters == rhs.m_oParameters)
            && (*m_pSurface == *rhs.m_pSurface)
            && (m_oIdentifier == rhs.m_oIdentifier));
  }

  /**
   * @brief inequality operator
   *
   * @return @c true if both objects are not equal, otherwise @c false
   *
   * @sa Measurement::operator==
   */
  bool
  operator!=(const Measurement<Identifier, params...>& rhs) const
  {
    return !(*this == rhs);
  }

  static Projection_t
  projector()
  {
    return ParSet_t::projector();
  }

  friend std::ostream&
  operator<<(std::ostream& out, const Measurement<Identifier, params...>& m)
  {
    m.print(out);
    return out;
  }

protected:
  virtual std::ostream&
  print(std::ostream& out) const
  {
    out << sizeof...(params) << "D measurement: ";
    int dummy[sizeof...(params)] = {(out << params << ", ", 0)...};
    out << std::endl;
    out << "measured values:" << std::endl;
    out << parameters() << std::endl;
    out << "covariance matrix:" << std::endl;
    out << covariance() << std::endl;
    out << "at " << (associatedSurface().isFree() ? "free" : "non-free")
        << " surface:" << std::endl;
    out << associatedSurface();

    return out;
  }

private:
  ParSet_t       m_oParameters;  ///< measured parameter set
  const Surface* m_pSurface;  ///< surface at which the measurement took place
  Identifier     m_oIdentifier;  ///< identifier for this measurement
};

/**
 * @brief general type for any possible Measurement
 */
template <typename Identifier>
using FittableMeasurement =
    typename detail::fittable_type_generator<Identifier>::type;

struct SurfaceGetter : public boost::static_visitor<const Surface&>
{
public:
  template <typename Meas_t>
  const Surface&
  operator()(const Meas_t& m) const
  {
    return m.associatedSurface();
  }
};

template <BOOST_VARIANT_ENUM_PARAMS(typename T)>
const Surface&
getSurface(const boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>& m)
{
  static const SurfaceGetter sg = SurfaceGetter();
  return boost::apply_visitor(sg, m);
}
}  // end of namespace Acts

#endif  // ACTS_MEASUREMENT_H

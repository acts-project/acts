// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <memory>
#include <ostream>
#include <type_traits>
#include <utility>
#include "Acts/EventData/ParameterSet.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/fittable_type_generator.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

// forward declarations
class Surface;

/// @brief base class for Measurements
///
/// This class describes the measurement of track parameters at a certain
/// Surface
/// in the TrackingGeometry.
///
/// @note Identifier must be copy-constructible, move-constructible,
/// copy-assignable and move-assignable.
///
/// @tparam identifier_t
///
/// @test The behavior of this class is tested in the following unit test:
///       - \link Acts::Test::BOOST_AUTO_TEST_CASE(measurement_initialization)
/// initialization\endlink
///
/// @tparam Identifier identification object for this measurement
/// @tparam params     parameter pack containing the measured parameters
template <typename identifier_t, ParID_t... params>
class Measurement
{
  // check type conditions
  static_assert(std::is_copy_constructible<identifier_t>::value,
                "'Identifier' must be copy-constructible");
  static_assert(std::is_move_constructible<identifier_t>::value,
                "'Identifier' must be move-constructible");
  static_assert(std::is_copy_assignable<identifier_t>::value,
                "'Identifier' must be copy-assignable");
  static_assert(std::is_move_assignable<identifier_t>::value,
                "'Identifier' must be move-assignable");

private:
  // private typedef's

  /// type of the underlying ParameterSet object
  using ParSet_t = ParameterSet<params...>;

public:
  /// type of the vector containing the parameter values
  using ParVector_t = typename ParSet_t::ParVector_t;
  /// type of the covariance matrix of the measurement
  using CovMatrix_t = typename ParSet_t::CovMatrix_t;
  /// matrix type for projecting full parameter vector onto local parameter
  /// space
  using Projection_t = typename ParSet_t::Projection_t;

  /// @brief standard constructor
  ///
  /// Interface class for all possible measurements.
  ///
  /// @note Only a reference to the given surface is stored. The user must
  /// ensure
  /// that the lifetime of the @c Surface
  ///       object surpasses the lifetime of this Measurement object.<br />
  ///       The given parameter values are interpreted as values to the
  /// parameters as defined in the class template
  ///       argument @c params.
  ///
  /// @attention The current design will fail if the in-memory location of
  /// the @c
  /// Surface object is changed (e.g.
  ///            if it is stored in a container and this gets relocated).
  ///
  /// @param surface surface at which the measurement took place
  /// @param id identification object for this measurement
  /// @param cov covariance matrix of the measurement.
  /// @param head,values consistent number of parameter values of the
  /// measurement
  template <typename... Tail>
  Measurement(std::shared_ptr<const Surface> surface,
              const identifier_t&            id,
              CovMatrix_t                    cov,
              typename std::enable_if<sizeof...(Tail) + 1 == sizeof...(params),
                                      ParValue_t>::type head,
              Tail... values)
    : m_oParameters(std::make_unique<const CovMatrix_t>(std::move(cov)),
                    head,
                    values...)
    , m_pSurface(std::move(surface))
    , m_oIdentifier(id)
  {
    assert(m_pSurface);
  }

  /// @brief virtual destructor
  virtual ~Measurement() = default;

  /// @brief copy constructor
  ///
  /// @tparam identifier_t The identifier type
  /// @tparam params...The local parameter pack
  ///
  /// @param copy is the source for the copy
  Measurement(const Measurement<identifier_t, params...>& copy)
    : m_oParameters(copy.m_oParameters)
    , m_pSurface(copy.m_pSurface)
    , m_oIdentifier(copy.m_oIdentifier)
  {
  }

  /// @brief move constructor
  ///
  /// @tparam identifier_t The identifier type
  /// @tparam params...The local parameter pack
  ///
  /// @param rhs is the source for the move
  Measurement(Measurement<identifier_t, params...>&& rhs)
    : m_oParameters(std::move(rhs.m_oParameters))
    , m_pSurface(std::move(rhs.m_pSurface))
    , m_oIdentifier(std::move(rhs.m_oIdentifier))
  {
  }

  /// @brief copy assignment operator
  ///
  /// @tparam identifier_t The identifier type
  /// @tparam params...The local parameter pack
  ///
  /// @param rhs is the source for the assignment
  Measurement<identifier_t, params...>&
  operator=(const Measurement<identifier_t, params...>& rhs)
  {
    m_oParameters = rhs.m_oParameters;
    m_pSurface    = rhs.m_pSurface;
    m_oIdentifier = rhs.m_oIdentifier;

    return *this;
  }

  /// @brief move assignment operator
  ///
  /// @tparam identifier_t The identifier type
  /// @tparam params...The local parameter pack
  ///
  /// @param rhs is the source for the move assignment
  Measurement<identifier_t, params...>&
  operator=(Measurement<identifier_t, params...>&& rhs)
  {
    m_oParameters = std::move(rhs.m_oParameters);
    m_pSurface    = std::move(rhs.m_pSurface);
    m_oIdentifier = std::move(rhs.m_oIdentifier);

    return *this;
  }

  /// @brief retrieve stored value for given parameter
  ///
  /// @tparam parameter identifier for the parameter to be retrieved
  /// @remark @c parameter must be part of the template parameter pack @c
  /// params.
  /// Otherwise a compile-time
  ///         error is generated.
  ///
  /// @return value of the stored parameter
  template <ParID_t parameter>
  ParValue_t
  get() const
  {
    return m_oParameters.template getParameter<parameter>();
  }

  /// @brief access vector with measured parameter values
  ///
  /// @return column vector whose size is equal to the dimensionality of this
  /// Measurement. The values are
  ///         given for the measured parameters in the order defined by the
  ///         class
  /// template argument @c params.
  ParVector_t
  parameters() const
  {
    return m_oParameters.getParameters();
  }

  /// @brief access covariance matrix of the measured parameter values
  ///
  /// @return covariance matrix of the measurement
  CovMatrix_t
  covariance() const
  {
    return *m_oParameters.getCovariance();
  }

  /// @brief retrieve stored uncertainty for given parameter
  ///
  /// @tparam parameter identifier for the parameter to be retrieved
  /// @remark @c parameter must be part of the template parameter pack @c
  /// params.
  /// Otherwise a compile-time
  ///         error is generated.
  ///
  /// @return uncertainty \f$\sigma \ge 0\f$ for given parameter
  template <ParID_t parameter>
  ParValue_t
  uncertainty() const
  {
    return m_oParameters.template getUncertainty<parameter>();
  }

  /// @brief number of measured parameters
  ///
  /// @return number of measured parameters
  static constexpr unsigned int
  size()
  {
    return ParSet_t::size();
  }

  /// @brief access associated surface
  ///
  /// @pre The @c Surface object used to construct this @c Measurement object
  /// must still be valid
  ///      at the same memory location.
  ///
  /// @return reference to surface at which the measurement took place
  const Acts::Surface&
  referenceSurface() const
  {
    return *m_pSurface;
  }

  /// @brief access to global measurement identifier
  ///
  /// @return identifier object
  identifier_t
  identifier() const
  {
    return m_oIdentifier;
  }

  /// @brief calculate residual with respect to given track parameters
  ///
  /// @note It is checked that the residual for non-local parameters are in
  /// valid
  /// range (e.g.
  ///       residuals in \f$\phi\f$ are corrected).
  ///
  /// @todo Implement check that TrackParameters are defined at the same Surface
  /// as the Measurement is.
  /// @todo Implement validity check for residuals of local parameters.
  ///
  /// @param trackPars reference TrackParameters object
  ///
  /// @return vector with the residual parameter values (in valid range)
  ///
  /// @sa ParameterSet::residual
  ParVector_t
  residual(const TrackParameters& trackPars) const
  {
    return m_oParameters.residual(trackPars.getParameterSet());
  }

  /// @brief equality operator
  ///
  /// @return @c true if parameter sets and associated surfaces compare equal,
  /// otherwise @c false
  virtual bool
  operator==(const Measurement<identifier_t, params...>& rhs) const
  {
    return ((m_oParameters == rhs.m_oParameters)
            && (*m_pSurface == *rhs.m_pSurface)
            && (m_oIdentifier == rhs.m_oIdentifier));
  }

  /// @brief inequality operator
  ///
  /// @return @c true if both objects are not equal, otherwise @c false
  ///
  /// @sa Measurement::operator==
  bool
  operator!=(const Measurement<identifier_t, params...>& rhs) const
  {
    return !(*this == rhs);
  }

  /// @projection operator
  static Projection_t
  projector()
  {
    return ParSet_t::projector();
  }

  friend std::ostream&
  operator<<(std::ostream& out, const Measurement<identifier_t, params...>& m)
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
    dummy[0]                     = dummy[0];
    out << std::endl;
    out << "measured values:" << std::endl;
    out << parameters() << std::endl;
    out << "covariance matrix:" << std::endl;
    out << covariance() << std::endl;
    out << "at " << (referenceSurface().isFree() ? "free" : "non-free")
        << " surface:" << std::endl;
    out << referenceSurface();

    return out;
  }

private:
  ParSet_t m_oParameters;  ///< measured parameter set
  std::shared_ptr<const Surface>
               m_pSurface;     ///< surface at which the measurement took place
  identifier_t m_oIdentifier;  ///< identifier for this measurement
};

/// @brief FittableMeasurement boost_variant type
template <typename identifier_t>
using FittableMeasurement =
    typename detail::fittable_type_generator<identifier_t>::type;

}  // namespace Acts

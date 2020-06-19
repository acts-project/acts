// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
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
#include "Acts/EventData/SourceLinkConcept.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/fittable_type_generator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

// forward declarations
class Surface;

/// @brief base class for Measurements
///
/// This class describes the measurement of track parameters at a certain
/// Surface in the TrackingGeometry.
///
/// The measurement is in local parameters and will not provide localToGlobal
/// information. It is thus free from any Context.
///
/// @tparam source_link_t the templated class that allows to link back to
/// the source used to create this measurement, this can simply be an identifier
/// or an extended object that allows to navigate back to the original source.
/// The source link is necessary when e.g. calibrating the object or in other
/// circumstances where the pure mathematical view of a measurement is not
/// sufficient.
/// The name _link_ indicates that there is a relationship to a unique source
/// object, i.e. if a measurement appears on several tracks, then there must
/// be a uniquely identifyable common source for those.
///
/// @test The behavior of this class is tested in the following unit test:
///       - \link Acts::Test::BOOST_AUTO_TEST_CASE(measurement_initialization)
/// initialization\endlink
///
/// @tparam Identifier identification object for this measurement
/// @tparam params     parameter pack containing the measured parameters
template <typename source_link_t, ParID_t... params>
class Measurement {
  // check type conditions
  static_assert(SourceLinkConcept<source_link_t>,
                "Source link does not fulfill SourceLinkConcept");

 private:
  // private typedefs

  /// type of the underlying ParameterSet object
  using ParSet_t = ParameterSet<BoundParametersIndices, params...>;

 public:
  /// type of the vector containing the parameter values
  using ParameterVector = typename ParSet_t::ParameterVector;
  /// type of the covariance matrix of the measurement
  using CovarianceMatrix = typename ParSet_t::CovarianceMatrix;
  /// matrix type for projecting full parameter vector onto local parameters
  using Projection = typename ParSet_t::Projection;

  /// Delete the default constructor
  Measurement() = delete;

  /// @brief standard constructor
  ///
  /// Concrete class for all possible measurements.
  ///
  /// @note Only a reference to the given surface is stored. The user must
  /// ensure that the lifetime of the @c Surface object surpasses the lifetime
  /// of this Measurement object.
  /// The given parameter values are interpreted as values to the
  /// parameters as defined in the class template argument @c params.
  ///
  /// @attention The current design will fail if the in-memory location of
  /// the @c Surface object is changed (e.g. if it is stored in a
  /// container and this gets relocated).
  ///
  /// @param surface surface at which the measurement took place
  /// @param source object for this measurement
  /// @param cov covariance matrix of the measurement.
  /// @param head,values consistent number of parameter values of the
  /// measurement
  template <typename... Tail>
  Measurement(std::shared_ptr<const Surface> surface,
              const source_link_t& source, CovarianceMatrix cov,
              typename std::enable_if<sizeof...(Tail) + 1 == sizeof...(params),
                                      ParValue_t>::type head,
              Tail... values)
      : m_oParameters(std::move(cov), head, values...),
        m_pSurface(std::move(surface)),
        m_sourceLink(source) {
    assert(m_pSurface);
  }

  /// @brief virtual destructor
  virtual ~Measurement() = default;

  /// @brief copy constructor
  ///
  /// @tparam source_link_t The identifier type
  /// @tparam params...The local parameter pack
  ///
  /// @param copy is the source for the copy
  Measurement(const Measurement<source_link_t, params...>& copy)
      : m_oParameters(copy.m_oParameters),
        m_pSurface(copy.m_pSurface),
        m_sourceLink(copy.m_sourceLink) {}

  /// @brief move constructor
  ///
  /// @tparam source_link_t The identifier type
  /// @tparam params...The local parameter pack
  ///
  /// @param other is the source for the move
  Measurement(Measurement<source_link_t, params...>&& other)
      : m_oParameters(std::move(other.m_oParameters)),
        m_pSurface(std::move(other.m_pSurface)),
        m_sourceLink(std::move(other.m_sourceLink)) {}

  /// @brief copy assignment operator
  ///
  /// @tparam source_link_t The identifier type
  /// @tparam params...The local parameter pack
  ///
  /// @param rhs is the source for the assignment
  Measurement<source_link_t, params...>& operator=(
      const Measurement<source_link_t, params...>& rhs) {
    // check for self-assignment
    if (&rhs != this) {
      m_oParameters = rhs.m_oParameters;
      m_pSurface = rhs.m_pSurface;
      m_sourceLink = rhs.m_sourceLink;
    }
    return *this;
  }

  /// @brief move assignment operator
  ///
  /// @tparam source_link_t The identifier type
  /// @tparam params...The local parameter pack
  ///
  /// @param rhs is the source for the move assignment
  Measurement<source_link_t, params...>& operator=(
      Measurement<source_link_t, params...>&& rhs) {
    m_oParameters = std::move(rhs.m_oParameters);
    m_pSurface = std::move(rhs.m_pSurface);
    m_sourceLink = std::move(rhs.m_sourceLink);
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
  ParValue_t get() const {
    return m_oParameters.template getParameter<parameter>();
  }

  /// @brief access vector with measured parameter values
  ///
  /// @return column vector whose size is equal to the dimensionality of this
  /// Measurement. The values are
  ///         given for the measured parameters in the order defined by the
  ///         class
  /// template argument @c params.
  ParameterVector parameters() const { return m_oParameters.getParameters(); }

  /// @brief access covariance matrix of the measured parameter values
  ///
  /// @return covariance matrix of the measurement
  CovarianceMatrix covariance() const { return *m_oParameters.getCovariance(); }

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
  ParValue_t uncertainty() const {
    return m_oParameters.template getUncertainty<parameter>();
  }

  /// @brief number of measured parameters
  ///
  /// @return number of measured parameters
  static constexpr unsigned int size() { return ParSet_t::size(); }

  /// @brief access associated surface
  ///
  /// @pre The @c Surface object used to construct this @c Measurement object
  /// must still be valid at the same memory location.
  ///
  /// @return reference to surface at which the measurement took place
  const Acts::Surface& referenceSurface() const { return *m_pSurface; }

  /// @brief link access to the source of the measurement.
  ///
  /// The source link can be simply an identifier or a more complicated
  /// object, see description above.
  ///
  /// @return source_link_t object
  const source_link_t& sourceLink() const { return m_sourceLink; }

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
  ParameterVector residual(const BoundParameters& trackPars) const {
    return m_oParameters.residual(trackPars.getParameterSet());
  }

  /// @brief equality operator
  ///
  /// @return @c true if parameter sets and associated surfaces compare equal,
  /// otherwise @c false
  virtual bool operator==(
      const Measurement<source_link_t, params...>& rhs) const {
    return ((m_oParameters == rhs.m_oParameters) &&
            (*m_pSurface == *rhs.m_pSurface) &&
            (m_sourceLink == rhs.m_sourceLink));
  }

  /// @brief inequality operator
  ///
  /// @return @c true if both objects are not equal, otherwise @c false
  ///
  /// @sa Measurement::operator==
  bool operator!=(const Measurement<source_link_t, params...>& rhs) const {
    return !(*this == rhs);
  }

  /// @projection operator
  static Projection projector() { return ParSet_t::projector(); }

  friend std::ostream& operator<<(
      std::ostream& out, const Measurement<source_link_t, params...>& m) {
    m.print(out);
    return out;
  }

 protected:
  virtual std::ostream& print(std::ostream& out) const {
    out << sizeof...(params) << "D measurement: ";
    int dummy[sizeof...(params)] = {(out << params << ", ", 0)...};
    dummy[0] = dummy[0];
    out << std::endl;
    out << "measured values:" << std::endl;
    out << parameters() << std::endl;
    out << "covariance matrix:" << std::endl;
    out << covariance() << std::endl;
    return out;
  }

 private:
  ParSet_t m_oParameters;  ///< measured parameter set
  std::shared_ptr<const Surface>
      m_pSurface;  ///< surface at which the measurement took place

  source_link_t m_sourceLink;  ///< link to the source for this measurement
};

/**
 * Required factory metafunction which produces measurements.
 * This encodes the source_link_t and hides it from the type generator.
 */
template <typename source_link_t>
struct fittable_measurement_helper {
  template <Acts::ParID_t... pars>
  struct meas_factory {
    using type = Measurement<source_link_t, pars...>;
  };

  using type =
      typename detail::type_generator_t<meas_factory, eBoundParametersSize>;
};

/**
 * @brief FittableMeasurement variant type
 */
template <typename source_link_t>
using FittableMeasurement =
    typename fittable_measurement_helper<source_link_t>::type;

}  // namespace Acts

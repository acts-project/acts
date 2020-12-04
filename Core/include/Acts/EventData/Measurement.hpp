// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/ParameterSet.hpp"
#include "Acts/EventData/SourceLinkConcept.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/CalculateResiduals.hpp"
#include "Acts/EventData/detail/ParameterTraits.hpp"
#include "Acts/EventData/detail/PrintParameters.hpp"
#include "Acts/EventData/detail/fittable_type_generator.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/detail/Subspace.hpp"

#include <array>
#include <iosfwd>
#include <memory>
#include <ostream>
#include <type_traits>
#include <utility>
#include <variant>

namespace Acts {

// forward declarations
class Surface;

namespace detail {
/// @brief Deduction of the measuring geometry object based on the used indices
template <typename T>
struct ReferenceObject {};
template <>
struct ReferenceObject<BoundIndices> {
  using type = Surface;
};
template <>
struct ReferenceObject<FreeIndices> {
  using type = Volume;
};
}  // namespace detail

/// @brief base class for Measurements
///
/// This class describes the measurement of track parameters at a certain
/// Surface/Volume in the TrackingGeometry.
///
/// If the measurement is in local parameters and will not provide localToGlobal
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
/// @tparam parameter_indices_t Enum of parameter identifier
/// @tparam params     parameter pack containing the measured parameters
template <typename source_link_t, typename parameter_indices_t,
          parameter_indices_t... params>
class Measurement {
  // check type conditions
  static_assert(SourceLinkConcept<source_link_t>,
                "Source link does not fulfill SourceLinkConcept");

  // type of the underlying ParameterSet object
  using ParamSet = ParameterSet<parameter_indices_t, params...>;

 public:
  using Scalar = typename ParamSet::Scalar;
  /// type of the vector containing the parameter values
  using ParametersVector = typename ParamSet::ParametersVector;
  /// Vector type containing all parameters from the same space
  using FullParametersVector = typename ParamSet::FullParametersVector;
  /// type of the covariance matrix of the measurement
  using CovarianceMatrix = typename ParamSet::CovarianceMatrix;
  /// matrix type for projecting full parameter vector onto measured parameters
  using ProjectionMatrix = typename ParamSet::ProjectionMatrix;
  /// Object type that corresponds to the measurement
  using RefObject = typename detail::ReferenceObject<parameter_indices_t>::type;

  /// Delete the default constructor
  Measurement() = delete;

  /// @brief standard constructor for surface/volume measurements
  ///
  /// Concrete class for all possible measurements.
  ///
  /// @param referenceObject surface/volume origin of the measurement
  /// @param source object for this measurement
  /// @param cov covariance matrix of the measurement.
  /// @param head,values consistent number of parameter values of the
  /// measurement
  template <typename... Tail>
  Measurement(std::shared_ptr<const RefObject> referenceObject,
              const source_link_t& source, CovarianceMatrix cov,
              typename std::enable_if<sizeof...(Tail) + 1 == sizeof...(params),
                                      Scalar>::type head,
              Tail... values)
      : m_oParameters(std::move(cov), head, values...),
        m_pReferenceObject(std::move(referenceObject)),
        m_sourceLink(source) {
    assert(m_pReferenceObject);
  }

  /// @brief standard constructor for surface/volume measurements
  ///
  /// Concrete class for all possible measurements, built from properly
  /// formatted covariance matrix and parameter vector
  ///
  /// @param referenceObject surface/volume origin of the measurement
  /// @param source object for this measurement
  /// @param cov covariance matrix of the measurement
  /// @param vec parameter vector of the measurement
  Measurement(std::shared_ptr<const RefObject> referenceObject,
              const source_link_t& source, CovarianceMatrix cov,
              ParametersVector vec)
      : m_oParameters(std::move(cov), std::move(vec)),
        m_pReferenceObject(std::move(referenceObject)),
        m_sourceLink(source) {
    assert(m_pReferenceObject);
  }

  /// @brief standard constructor for surface/volume measurements
  ///
  /// Concrete class for all possible measurements, built from properly
  /// formatted parameter set
  ///
  /// @param referenceObject surface/volume origin of the measurement
  /// @param source object for this measurement
  /// @param pset parameter set of the measurement
  Measurement(std::shared_ptr<const RefObject> referenceObject,
              const source_link_t& source, ParamSet pset)
      : m_oParameters(std::move(pset)),
        m_pReferenceObject(std::move(referenceObject)),
        m_sourceLink(source) {
    assert(m_pReferenceObject);
  }

  virtual ~Measurement() = default;

  /// @brief copy constructor
  ///
  /// @tparam source_link_t The identifier type
  /// @tparam parameter_indices_t Enum of parameter identifier
  /// @tparam params...The local parameter pack
  ///
  /// @param copy is the source for the copy
  Measurement(
      const Measurement<source_link_t, parameter_indices_t, params...>& copy)
      : m_oParameters(copy.m_oParameters),
        m_pReferenceObject(copy.m_pReferenceObject),
        m_sourceLink(copy.m_sourceLink) {}

  /// @brief move constructor
  ///
  /// @tparam source_link_t The identifier type
  /// @tparam parameter_indices_t Enum of parameter identifier
  /// @tparam params...The local parameter pack
  ///
  /// @param other is the source for the move
  Measurement(
      Measurement<source_link_t, parameter_indices_t, params...>&& other)
      : m_oParameters(std::move(other.m_oParameters)),
        m_pReferenceObject(std::move(other.m_pReferenceObject)),
        m_sourceLink(std::move(other.m_sourceLink)) {}

  /// @brief copy assignment operator
  ///
  /// @tparam source_link_t The identifier type
  /// @tparam parameter_indices_t Enum of parameter identifier
  /// @tparam params...The local parameter pack
  ///
  /// @param rhs is the source for the assignment
  Measurement<source_link_t, parameter_indices_t, params...>& operator=(
      const Measurement<source_link_t, parameter_indices_t, params...>& rhs) {
    // check for self-assignment
    if (&rhs != this) {
      m_oParameters = rhs.m_oParameters;
      m_pReferenceObject = rhs.m_pReferenceObject;
      m_sourceLink = rhs.m_sourceLink;
    }
    return *this;
  }

  /// @brief move assignment operator
  ///
  /// @tparam source_link_t The identifier type
  /// @tparam parameter_indices_t Enum of parameter identifier
  /// @tparam params...The local parameter pack
  ///
  /// @param rhs is the source for the move assignment
  Measurement<source_link_t, parameter_indices_t, params...>& operator=(
      Measurement<source_link_t, parameter_indices_t, params...>&& rhs) {
    m_oParameters = std::move(rhs.m_oParameters);
    m_pReferenceObject = std::move(rhs.m_pReferenceObject);
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
  template <parameter_indices_t parameter>
  Scalar get() const {
    return m_oParameters.template getParameter<parameter>();
  }

  /// @brief checks whether a given parameter is included in this set of
  ///  parameters
  ///
  /// @tparam parameter identifier for the parameter to be retrieved
  /// @remark @c parameter must be part of the template parameter pack @c
  /// params. Otherwise a compile-time error is generated.
  ///
  /// @return @c true if the parameter is stored in this set, otherwise @c false
  template <parameter_indices_t parameter>
  bool contains() const {
    return m_oParameters.template contains<parameter>();
  }

  /// @brief access vector with measured parameter values
  ///
  /// @return column vector whose size is equal to the dimensionality of this
  /// Measurement. The values are given for the measured parameters in the
  /// order defined by the class template argument @c params.
  ///
  /// @return A pure vector type of length size_of(params...)
  const ParametersVector& parameters() const {
    return m_oParameters.getParameters();
  }

  /// @brief access covariance matrix of the measured parameter values
  ///
  /// @return covariance matrix of the measurement
  const CovarianceMatrix& covariance() const {
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
  template <parameter_indices_t parameter>
  Scalar uncertainty() const {
    return m_oParameters.template getUncertainty<parameter>();
  }

  /// @brief number of measured parameters
  ///
  /// @return number of measured parameters
  static constexpr unsigned int size() { return ParamSet::size(); }

  /// @brief access associated object
  ///
  /// @pre The @c ReferenceObject object used to construct this @c Measurement
  /// object must still be valid at the same memory location.
  ///
  /// @return reference to surface/volume associated which the measurement
  const RefObject& referenceObject() const { return *m_pReferenceObject; }

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
  /// @todo Implement check that TrackParameters are defined at the same
  /// Surface/Volume as the Measurement is.
  /// @todo Implement validity check for residuals of parameters.
  ///
  /// @param trackPars reference TrackParameters object
  ///
  /// @return vector with the residual parameter values (in valid range)
  ///
  /// @sa ParameterSet::residual
  ParametersVector residual(const FullParametersVector& trackPars) const {
    return m_oParameters.residual(trackPars);
  }

  /// @brief equality operator
  ///
  /// @return @c true if parameter sets and associated surfaces/volumes compare
  /// equal, otherwise @c false
  virtual bool operator==(const Measurement<source_link_t, parameter_indices_t,
                                            params...>& rhs) const {
    return ((m_oParameters == rhs.m_oParameters) &&
            (m_pReferenceObject == rhs.m_pReferenceObject) &&
            (m_sourceLink == rhs.m_sourceLink));
  }

  /// @brief inequality operator
  ///
  /// @return @c true if both objects are not equal, otherwise @c false
  ///
  /// @sa Measurement::operator==
  bool operator!=(const Measurement<source_link_t, parameter_indices_t,
                                    params...>& rhs) const {
    return !(*this == rhs);
  }

  /// @ProjectionMatrix operator
  static const ProjectionMatrix& projector() { return ParamSet::projector(); }

  friend std::ostream& operator<<(
      std::ostream& out,
      const Measurement<source_link_t, parameter_indices_t, params...>& m) {
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
  ParamSet m_oParameters;  ///< measured parameter set
  std::shared_ptr<const RefObject> m_pReferenceObject =
      nullptr;                 ///< object which corresponds to the measurement
  source_link_t m_sourceLink;  ///< link to the source for this measurement
};

/**
 * Required factory metafunction which produces measurements.
 * This encodes the source_link_t and hides it from the type generator.
 */
template <typename source_link_t>
struct fittable_measurement_helper {
  template <BoundIndices... pars>
  struct meas_factory {
    using type = Measurement<source_link_t, BoundIndices, pars...>;
  };

  using type = typename detail::type_generator_t<BoundIndices, meas_factory>;
};

template <typename source_link_t>
struct fittable_volume_measurement_helper {
  template <FreeIndices... pars>
  struct meas_factory {
    using type = Measurement<source_link_t, FreeIndices, pars...>;
  };

  using type = typename detail::type_generator_t<FreeIndices, meas_factory>;
};

/**
 * @brief Measurement variant types
 */
template <typename source_link_t>
using FittableMeasurement =
    typename fittable_measurement_helper<source_link_t>::type;

template <typename source_link_t>
using FittableVolumeMeasurement =
    typename fittable_volume_measurement_helper<source_link_t>::type;

/// A measurement of a fixed-size subspace of the full parameters.
///
/// @tparam source_link_t Source link type to connect to the detector readout
/// @tparam indices_t Parameter index type, determines the full parameter space
/// @tparam kSize Size of the parameter subset.
///
/// The measurement is templated on the measurement subspace size and its
/// content (indices) is defined at runtime. In previous implementations,
/// the measurement was templated directly on the subspace indices. This lead to
/// a combinatorial explosion and severe penalty during compilation since all
/// possible combinations need to be considered for the measurement variant.
/// However, for the computations only the size of the measurement subspace is
/// relevant.
///
/// The measurement intenionally does not store a pointer/reference to the
/// explicit reference object in the geometry hierarchy, i.e. the surface or
/// volume. The reference object can already be identified via the geometry
/// identifier provided by the source link. Since a measurement
/// **must** be anchored within the geometry hierarchy, all measurement surfaces
/// and volumes **must** provide valid geometry identifiers. In all use-cases,
/// e.g. Kalman filtering, a pointer/reference to the reference object is
/// available before the measurement is accessed; e.g. the propagator provides
/// the surface pointer during navigation, which is then used to lookup possible
/// measurements.
///
/// The pointed-to geometry object would differ depending on the parameter type.
/// This means either, that there needs to be an additional variable type or
/// that a pointer to a base object is stored (requiring a `dynamic_cast` later
/// on). Both variants add additional complications. Since the geometry object
/// is not required anyways (as discussed above), not storing it removes all
/// these complications altogether.
template <typename source_link_t, typename indices_t, size_t kSize>
class FixedSizeMeasurement {
  static_assert(SourceLinkConcept<source_link_t>,
                "Source link does not fulfill SourceLinkConcept");

  static constexpr size_t kFullSize = detail::kParametersSize<indices_t>;

  using Subspace = detail::FixedSizeSubspace<kFullSize, kSize>;

 public:
  using Scalar = typename detail::ParametersScalar<indices_t>;
  /// Vector type containing for measured parameter values.
  using ParametersVector = ActsVector<Scalar, kSize>;
  /// Matrix type for the measurement covariance.
  using CovarianceMatrix = ActsSymMatrix<Scalar, kSize>;
  /// Vector type containing all parameters in the same space.
  using FullParametersVector = ActsVector<Scalar, kFullSize>;
  using ProjectionMatrix = ActsMatrix<Scalar, kSize, kFullSize>;
  using ExpansionMatrix = ActsMatrix<Scalar, kFullSize, kSize>;

  /// Construct from source link, subset indices, and measured data.
  ///
  /// @tparam index_container_t Forward-iterable Container for indices
  /// @tparam parameters_t Input parameters vector type
  /// @tparam covariance_t Input covariance matrix type
  /// @param source The link that connects to the underlying detector readout
  /// @param indices Which parameters are measured
  /// @param params Measured parameters values
  /// @param cov Measured parameters covariance
  ///
  /// @note The indices must be ordered and must describe/match the content
  ///   of parameters and covariance.
  template <typename index_container_t, typename parameters_t,
            typename covariance_t>
  FixedSizeMeasurement(source_link_t source, const index_container_t& indices,
                       const Eigen::MatrixBase<parameters_t>& params,
                       const Eigen::MatrixBase<covariance_t>& cov)
      : m_source(std::move(source)),
        m_subspace(indices),
        m_params(params),
        m_cov(cov) {
    // TODO we should be able to support arbitrary ordering, by sorting the
    //   indices and reordering parameters/covariance. since the parameter order
    //   can be modified by the user, the user can not always know what the
    //   right order is. another option is too remove the requirement for index
    //   ordering from the subspace types, but that will make it harder to
    //   refactor their implementation later on.
  }
  /// A measurement can only be constructed with valid parameters.
  FixedSizeMeasurement() = delete;
  FixedSizeMeasurement(const FixedSizeMeasurement&) = default;
  FixedSizeMeasurement(FixedSizeMeasurement&&) = default;
  ~FixedSizeMeasurement() = default;
  FixedSizeMeasurement& operator=(const FixedSizeMeasurement&) = default;
  FixedSizeMeasurement& operator=(FixedSizeMeasurement&&) = default;

  /// Source link that connects to the underlying detector readout.
  const source_link_t& sourceLink() const { return m_source; }

  /// Number of measured parameters.
  static constexpr size_t size() { return kSize; }

  /// Check if a specific parameter is part of this measurement.
  bool contains(indices_t i) const { return m_subspace.contains(i); }

  /// Measured parameters values.
  const ParametersVector& parameters() const { return m_params; }

  /// Measured parameters covariance.
  const CovarianceMatrix& covariance() const { return m_cov; }

  /// Projection matrix from the full space into the measured subspace.
  ProjectionMatrix projector() const {
    return m_subspace.template projector<Scalar>();
  }

  /// Expansion matrix from the measured subspace into the full space.
  ///
  /// This is equivalent to the tranpose of the projection matrix only in the
  /// case of a trivial projection matrix. While this is the case here, it is
  /// still recommended to use the expansion matrix directly in cases where it
  /// is explicitely used.
  ExpansionMatrix expander() const {
    return m_subspace.template expander<Scalar>();
  }

  /// Compute residuals in the measured subspace.
  ///
  /// @param reference Reference parameters in the full space.
  ///
  /// This computes the difference `measured - reference` taking into account
  /// the allowed parameter ranges. Only the reference values in the measured
  /// subspace are used for the computation.
  ParametersVector residuals(const FullParametersVector& reference) const {
    ParametersVector res = ParametersVector::Zero();
    detail::calculateResiduals(static_cast<indices_t>(kSize),
                               m_subspace.indices(), reference, m_params, res);
    return res;
  }

  std::ostream& operator<<(std::ostream& os) const {
    detail::printMeasurement(os, static_cast<indices_t>(kSize),
                             m_subspace.indices().data(), m_params, m_cov);
    return os;
  }

 private:
  source_link_t m_source;
  Subspace m_subspace;
  ParametersVector m_params;
  CovarianceMatrix m_cov;
};

/// Construct a fixed-size measurement for the given indices.
///
/// @tparam source_link_t Source link type to connect to the detector readout
/// @tparam parameters_t Input parameters vector type
/// @tparam covariance_t Input covariance matrix type
/// @tparam indices_t Parameter index type, determines the full parameter space
/// @tparam tail_indices_t Helper types required to support variadic arguments;
///   all types must be convertibale to `indices_t`.
/// @param source The link that connects to the underlying detector readout
/// @param params Measured parameters values
/// @param cov Measured parameters covariance
/// @param index0 Required parameter index, measurement must be at least 1d
/// @param tailIndices Additional parameter indices for larger measurements
/// @return Fixed-size measurement w/ the correct type and given inputs
///
/// This helper function can be used to create a fixed-size measurement using an
/// explicit set of indices, e.g.
///
///     auto m = makeMeasurement(s, p, c, eBoundLoc0, eBoundTime);
///
/// for a 2d measurement w/ one position and time.
///
/// @note The indices must be ordered and must be consistent with the content of
///   parameters and covariance.
template <typename source_link_t, typename parameters_t, typename covariance_t,
          typename indices_t, typename... tail_indices_t>
auto makeMeasurement(source_link_t source,
                     const Eigen::MatrixBase<parameters_t>& params,
                     const Eigen::MatrixBase<covariance_t>& cov,
                     indices_t index0, tail_indices_t... tailIndices)
    -> FixedSizeMeasurement<source_link_t, indices_t,
                            1u + sizeof...(tail_indices_t)> {
  using IndexContainer = std::array<indices_t, 1u + sizeof...(tail_indices_t)>;
  return {std::move(source), IndexContainer{index0, tailIndices...}, params,
          cov};
}

namespace detail {
// Recursive construction of the measurement variant. `kN` is counted down until
// zero while the sizes are accumulated in the parameter pack.
//
// Example:
//
//        VariantMeasurementGenerator<..., 4>
//     -> VariantMeasurementGenerator<..., 3, 4>
//     -> VariantMeasurementGenerator<..., 2, 3, 4>
//     -> VariantMeasurementGenerator<..., 1, 2, 3, 4>
//     -> VariantMeasurementGenerator<..., 0, 1, 2, 3, 4>
//
template <typename source_link_t, typename indices_t, size_t kN,
          size_t... kSizes>
struct VariantMeasurementGenerator
    : VariantMeasurementGenerator<source_link_t, indices_t, kN - 1u, kN,
                                  kSizes...> {};
template <typename source_link_t, typename indices_t, size_t... kSizes>
struct VariantMeasurementGenerator<source_link_t, indices_t, 0u, kSizes...> {
  using Type =
      std::variant<FixedSizeMeasurement<source_link_t, indices_t, kSizes>...>;
};
}  // namespace detail

/// Variant that can contain all possible measurements in a parameter space.
///
/// @tparam source_link_t Source link type to connect to the detector readout
/// @tparam indices_t Parameter index type, determines the full parameter space
template <typename source_link_t, typename indices_t>
using VariantMeasurement = typename detail::VariantMeasurementGenerator<
    source_link_t, indices_t, detail::kParametersSize<indices_t>>::Type;

/// Variant that can hold all possible bound measurements.
///
/// @tparam source_link_t Source link type to connect to the detector readout
template <typename source_link_t>
using BoundVariantMeasurement = VariantMeasurement<source_link_t, BoundIndices>;

/// Variant that can hold all possible free measurements.
///
/// @tparam source_link_t Source link type to connect to the detector readout
template <typename source_link_t>
using FreeVariantMeasurement = VariantMeasurement<source_link_t, FreeIndices>;

template <typename source_link_t, typename indices_t>
std::ostream& operator<<(
    std::ostream& os, const VariantMeasurement<source_link_t, indices_t>& vm) {
  return std::visit([&](const auto& m) { return (os << m); }, vm);
}

}  // namespace Acts

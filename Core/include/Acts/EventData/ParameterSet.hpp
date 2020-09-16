// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/detail/DifferenceCalculator.hpp"
#include "Acts/EventData/detail/ValueCorrector.hpp"
#include "Acts/EventData/detail/full_parameter_set.hpp"
#include "Acts/EventData/detail/initialize_parameter_set.hpp"
#include "Acts/EventData/detail/make_projection_matrix.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/detail/MPL/are_sorted.hpp"
#include "Acts/Utilities/detail/MPL/are_within.hpp"
#include "Acts/Utilities/detail/MPL/at_index.hpp"
#include "Acts/Utilities/detail/MPL/get_position.hpp"
#include "Acts/Utilities/detail/MPL/is_contained.hpp"

#include <optional>
#include <type_traits>

namespace Acts {

// Parameter sets corresponding to the full parameters vector
using FullBoundParameterSet = typename detail::full_parset<BoundIndices>::type;
using FullFreeParameterSet = typename detail::full_parset<FreeIndices>::type;

/**
 * @class ParameterSet
 *
 * @brief Description of a set of (local) parameters
 *
 * @pre
 * The template parameter @c ParameterPolicy must fulfill the following
 * requirements:
 *  -# It must contain a <tt>typedef #parameter_indices_t</tt> specifying an
 * integral type used to identify different parameters. This could for example
 * be an @c enum, @c short, or <tt>unsigned int</tt>. This @c typedef must be
 * convertible to an <tt>unsigned int</tt>
 *  -# It must contain a <tt>typedef #Scalar</tt> specifying the type of
 * the parameter values. This could for
 *     instance be @c double, or @c float.
 *  -# It must contain a definition of an integral constant named @c N which is
 * assignable to an <tt>unsigned
 *     int</tt> and which is equal to the total number of parameters in the
 * system.
 * @pre
 *
 * The template parameter pack @c params must be given in a strictly ascending
 * order. The parameter pack must
 * be non-empty and it cannot contain more elements than
 * <tt>parsSize</tt>.
 *
 * @test The behavior of this class is tested in the following unit tests:
 *       - \link Acts::Test::BOOST_AUTO_TEST_CASE(parset_consistency_tests)
 * general consistency\endlink
 *       - \link Acts::Test::BOOST_AUTO_TEST_CASE(parset_copy_assignment_tests)
 * copy/assignment/swap\endlink
 *       - \link Acts::Test::BOOST_AUTO_TEST_CASE(parset_comparison_tests)
 * comparison operators\endlink
 *       - \link Acts::Test::BOOST_AUTO_TEST_CASE(parset_projection_tests)
 * projection matrices\endlink
 *       - \link Acts::Test::BOOST_AUTO_TEST_CASE(parset_residual_tests)
 * residual calculation\endlink
 *
 * @tparam ParameterPolicy  struct or class containing the parameter definitions
 * (see above)
 * @tparam params           parameter pack containing the (local) parameters
 * stored in this class
 */
template <typename parameter_indices_t, parameter_indices_t... params>
class ParameterSet {
 private:
  /// number of parameters stored in this class
  static constexpr unsigned int kNumberOfParameters = sizeof...(params);
  /// Highest index in used parameter indices
  static constexpr unsigned int kSizeMax =
      detail::kParametersSize<parameter_indices_t>;

  // static assert to check that the template parameters are consistent
  static_assert(
      detail::are_sorted<true, true, parameter_indices_t, params...>::value,
      "parameter identifiers are not sorted");
  static_assert(
      detail::are_within<unsigned int, 0, kSizeMax,
                         static_cast<unsigned int>(params)...>::value,
      "parameter identifiers must be greater or "
      "equal to zero and smaller than the total number of parameters");
  static_assert(kNumberOfParameters > 0,
                "number of stored parameters can not be zero");
  static_assert(
      kNumberOfParameters <= kSizeMax,
      "number of stored parameters can not exceed number of total parameters");

 public:
  using Scalar = detail::ParametersScalar<parameter_indices_t>;
  using ParametersVector = ActsVector<Scalar, kNumberOfParameters>;
  /// Vector type containing all parameters from the same space
  using FullParametersVector = ActsVector<Scalar, kSizeMax>;
  using CovarianceMatrix = ActsSymMatrix<Scalar, kNumberOfParameters>;
  /// Projection matrix to project full parameters into the configured space.
  using ProjectionMatrix = ActsMatrix<Scalar, kNumberOfParameters, kSizeMax>;

  /**
   * @brief initialize values of stored parameters and their covariance matrix
   *
   * @note  No validation of the given covariance matrix is performed (e.g. that
   * it is symmetric).
   *
   * @param cov unique pointer to covariance matrix (nullptr is accepted)
   * @param head value for first parameter
   * @param values values for the remaining stored parameters
   */
  template <typename... Tail>
  ParameterSet(
      std::optional<CovarianceMatrix> cov,
      std::enable_if_t<sizeof...(Tail) + 1 == kNumberOfParameters, Scalar> head,
      Tail... values)
      : m_vValues(kNumberOfParameters) {
    if (cov) {
      m_optCovariance = std::move(*cov);
    }
    detail::initialize_parset<parameter_indices_t, params...>::init_vals(
        *this, head, values...);
  }

  /**
   * @brief initialize parameter values from vector and set their covariance
   * matrix
   *
   * @note The values in the passed vector are interpreted as parameter values
   * in the order given
   *       by the class template @c params. No validation of the given
   * covariance matrix is performed.
   *
   * @param cov unique pointer to covariance matrix (nullptr is accepted)
   * @param values vector with parameter values
   */
  ParameterSet(std::optional<CovarianceMatrix> cov,
               const ParametersVector& values)
      : m_vValues(kNumberOfParameters) {
    if (cov) {
      m_optCovariance = std::move(*cov);
    }
    detail::initialize_parset<parameter_indices_t, params...>::init_vec(*this,
                                                                        values);
  }

  // this class does not have a custom default constructor and thus should not
  // provide any custom default cstors, dstor, or assignment. see ISOCPP C.20.

  /**
   * @brief return index of parameter identifier in parameter list
   *
   * @tparam parameter identifier for the parameter to be retrieved
   * @remark @c parameter must be part of the template parameter pack @c params.
   *         Otherwise a compile-time error is generated.
   *
   * @return position of parameter in variadic template parameter set @c params
   */
  template <parameter_indices_t parameter>
  static constexpr size_t getIndex() {
    return detail::get_position<parameter_indices_t, parameter,
                                params...>::value;
  }

  /**
   * @brief return parameter identifier for given index
   *
   * @tparam index position of parameter identifier in @c params
   * @remark @c index must be a positive number smaller than the size of the
   *         parameter pack @c params. Otherwise a compile-time error is
   *         generated.
   *
   * @return parameter identifier at position @c index in variadic template
   *         parameter set @c params
   */
  template <size_t index>
  static constexpr parameter_indices_t getParameterIndex() {
    return detail::at_index<parameter_indices_t, index, params...>::value;
  }

  /**
   * @brief retrieve stored value for given parameter
   *
   * @tparam parameter identifier for the parameter to be retrieved
   * @remark @c parameter must be part of the template parameter pack @c params.
   *         Otherwise a compile-time error is generated.
   *
   * @return value of the stored parameter
   */
  template <parameter_indices_t parameter>
  Scalar getParameter() const {
    return m_vValues(getIndex<parameter>());
  }

  /**
   * @brief access vector with stored parameters
   *
   * @return column vector with @c #kNumberOfParameters rows
   */
  const ParametersVector& getParameters() const { return m_vValues; }

  /**
   * @brief sets value for given parameter
   *
   * @tparam parameter identifier for the parameter to be stored
   * @remark @c parameter must be part of the template parameter pack @c params.
   * Otherwise a compile-time
   *         error is generated.
   *
   * @return previously stored value of this parameter
   */
  template <parameter_indices_t parameter>
  void setParameter(Scalar value) {
    m_vValues(getIndex<parameter>()) =
        detail::ParameterTraits<parameter_indices_t, parameter>::getValue(
            value);
  }

  /**
   * @brief sets values of stored parameters
   *
   * The values of the given vector are interpreted as parameter values in the
   * order
   * of the class template `params...`.
   *
   * @param values vector of length #kNumberOfParameters
   */
  void setParameters(const ParametersVector& values) {
    detail::initialize_parset<parameter_indices_t, params...>::init_vec(*this,
                                                                        values);
  }

  /**
   * @brief checks whether a given parameter is included in this set of
   parameters

   * @tparam parameter identifier for the parameter to be retrieved
   * @remark @c parameter must be part of the template parameter pack @c params.
   Otherwise a compile-time
   *         error is generated.
   *
   * @return @c true if the parameter is stored in this set, otherwise @c false
   */
  template <parameter_indices_t parameter>
  bool contains() const {
    return detail::is_contained<parameter_indices_t, parameter,
                                params...>::value;
  }

  /**
   * @brief access covariance matrix for stored parameters
   *
   * @note The ownership of the covariance matrix is @b not transferred with
   * this call.
   *
   * @return raw pointer to covariance matrix (can be a nullptr)
   */
  const std::optional<CovarianceMatrix>& getCovariance() const {
    return m_optCovariance;
  }

  /**
   * @brief access uncertainty for individual parameter
   *
   * @tparam parameter identifier for the parameter to be retrieved
   * @remark @c parameter must be part of the template parameter pack @c params.
   * Otherwise a compile-time
   *         error is generated.
   *
   * @return uncertainty \f$\sigma \ge 0\f$ of given parameter, a negative value
   * is returned if no
   *         covariance matrix is set
   */
  template <parameter_indices_t parameter>
  Scalar getUncertainty() const {
    if (m_optCovariance) {
      size_t index = getIndex<parameter>();
      return sqrt((*m_optCovariance)(index, index));
    } else {
      return -1;
    }
  }

  /**
   * @brief update covariance matrix
   *
   * @note No validation of the given covariance matrix is performed.
   *
   * @param cov unique pointer to new covariance matrix (nullptr is accepted)
   */
  void setCovariance(const CovarianceMatrix& cov) { m_optCovariance = cov; }

  /**
   * @brief calculate residual difference to full parameter vector
   *
   * Calculate the residual differences of the stored parameter values with
   * respect to the corresponding parameter values in the full parameter vector.
   * Hereby, the residual vector is defined as
   *
   * \f[
   * \vec{r} = \left( \begin{array}{c} r_{i_1} \\ \vdots \\ r_{i_m} \end{array}
   * \right)
   *  = \left( \begin{array}{c} v_{i_1} \\ \vdots \\ v_{i_m} \end{array} \right)
   * -  \mathrm{Proj} \left( \begin{array}{c} v^0_{1} \\ \vdots \\ v^0_{N}
   * \end{array} \right)
   *  = \vec{v} - \mathrm{Proj} \left( \vec{v}^0 \right)
   * \f]
   *
   * where \f$\mathrm{Proj}\f$ is the projection matrix, \f$\vec{v}\f$ is the
   * vector of parameter values of this ParameterSet object and \f$\vec{v}^0\f$
   * is the full parameter value vector.
   *
   * @param boundParams Vector of bound parameters
   * @note Constraint and cyclic parameter value ranges of @p
   * BoundTrackParameters are not tested.
   * @note It is not tested whether @p BoundTrackParameters is at the same
   * reference object
   *
   * @return vector containing the residual parameter values of this
   * ParameterSet object
   *         with respect to the given full parameter vector
   *
   * @sa ParameterSet::projector
   */
  ParametersVector residual(const FullParametersVector& other) const {
    return detail::DifferenceCalculator<parameter_indices_t, params...>::run(
        m_vValues, projector() * other);
  }

  /**
   * @brief calculate residual difference to other parameter vector
   *
   * Calculate the residual differences of the stored parameter values with
   * respect to the values of
   * another ParameterSet object containing the same set of parameters. Hereby,
   * the residual vector is
   * defined as
   *
   * \f[
   * \vec{r} = \left( \begin{array}{c} r_{i_1} \\ \vdots \\ r_{i_m} \end{array}
   * \right)
   *  = \left( \begin{array}{c} v_{i_1} \\ \vdots \\ v_{i_m} \end{array} \right)
   * -  \left( \begin{array}{c} v^0_{1} \\ \vdots \\ v^0_{N} \end{array} \right)
   *  = \vec{v} - \left( \vec{v}^0 \right)
   * \f]
   *
   * where \f$\vec{v}\f$ is the vector of parameter values of this ParameterSet
   * object and \f$\vec{v}^0\f$
   * is the parameter value vector of the other ParameterSet object.
   *
   * @note Constraint and cyclic parameter value ranges are taken into account
   * when calculating
   *       the residual values.
   *
   * @param otherParSet ParameterSet object with identical set of contained
   * parameters
   *
   * @return vector containing the residual parameter values of this
   * ParameterSet object
   *         with respect to the given other parameter set
   */
  ParametersVector residual(const ParameterSet& otherParSet) const {
    return detail::DifferenceCalculator<parameter_indices_t, params...>::run(
        m_vValues, otherParSet.m_vValues);
  }

  /**
   * @brief get projection matrix
   *
   * The projection matrix performs a mapping of the full parameter space onto
   * the sub-space
   * spanned by the parameters defined in this ParameterSet object.
   *
   * @return constant matrix with @c #kNumberOfParameters rows and @c
   * #kSizeMax columns
   */
  static const ProjectionMatrix& projector() { return sProjector; }

  /**
   * @brief number of stored parameters
   *
   * @return number of stored parameters
   */
  static constexpr unsigned int size() { return kNumberOfParameters; }

  /**
   * @brief correct given parameter values
   *
   * Check that the given values are within in a valid range for the
   * corresponding parameter. If not, an
   * in-place correction is applied. The values are interpreted as parameter
   * values in the same order as
   * specified in the class template @c params.
   *
   * @param values vector with parameter values to be checked and corrected if
   * necessary
   */
  static void correctValues(ParametersVector& values) {
    detail::ValueCorrector<parameter_indices_t, params...>::run(values);
  }

  /// Compare with another parameter set for equality.
  bool operator==(const ParameterSet& other) const {
    return (this == &other) or ((m_vValues == other.m_vValues) and
                                (m_optCovariance == other.m_optCovariance));
  }
  /// Compare with another parameter set for inequality.
  bool operator!=(const ParameterSet& other) const { return !(*this == other); }

 private:
  /// column vector containing values of local parameters.
  ParametersVector m_vValues{ParametersVector::Zero()};
  /// optional covariance matrix.
  std::optional<CovarianceMatrix> m_optCovariance{std::nullopt};

  /// matrix to project full parameter vector onto local parameter space.
  static const ProjectionMatrix sProjector;
};

template <typename parameter_indices_t, parameter_indices_t... kParameters>
const typename ParameterSet<parameter_indices_t,
                            kParameters...>::ProjectionMatrix
    ParameterSet<parameter_indices_t, kParameters...>::sProjector =
        detail::make_projection_matrix<
            kSizeMax, static_cast<unsigned int>(kParameters)...>::init();

}  // namespace Acts

// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// STL include(s)
#include <memory>
#include <optional>
#include <type_traits>
#include <utility>

// Acts includes
#include "Acts/EventData/detail/full_parameter_set.hpp"
#include "Acts/EventData/detail/initialize_parameter_set.hpp"
#include "Acts/EventData/detail/make_projection_matrix.hpp"
#include "Acts/EventData/detail/residual_calculator.hpp"
#include "Acts/EventData/detail/value_corrector.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/detail/MPL/are_sorted.hpp"
#include "Acts/Utilities/detail/MPL/are_within.hpp"
#include "Acts/Utilities/detail/MPL/at_index.hpp"
#include "Acts/Utilities/detail/MPL/get_position.hpp"
#include "Acts/Utilities/detail/MPL/is_contained.hpp"

namespace Acts {
/// @cond
// forward type declaration for full parameter set
using FullParameterSet =
    typename detail::full_parset<BoundParametersIndices>::type;
using FullFreeParameterSet =
    typename detail::full_parset<FreeParametersIndices>::type;
/// @endcond

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
 *  -# It must contain a <tt>typedef #ParValue_t</tt> specifying the type of
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
  // local typedefs and constants
  using Self = ParameterSet<parameter_indices_t,
                            params...>;  ///< type of this parameter set
  static constexpr unsigned int kNumberOfParameters =
      sizeof...(params);  ///< number of parameters stored in this class
  static constexpr unsigned int kSizeMax = detail::ParametersSize<
      parameter_indices_t>::size;  ///< Highest index in used parameter indices

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
  // public typedefs
  /// matrix type for projecting full parameter vector onto local parameter
  /// space
  using Projection = ActsMatrix<ParValue_t, kNumberOfParameters, kSizeMax>;
  /// vector type for stored parameters
  using ParameterVector = ActsVector<ParValue_t, kNumberOfParameters>;
  /// type of covariance matrix
  using CovarianceMatrix = ActsSymMatrix<ParValue_t, kNumberOfParameters>;

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
      std::enable_if_t<sizeof...(Tail) + 1 == kNumberOfParameters, ParValue_t>
          head,
      Tail... values)
      : m_vValues(kNumberOfParameters) {
    if (cov) {
      m_optCovariance = std::move(*cov);
    }
    detail::initialize_parset<parameter_indices_t, params...>::init(*this, head,
                                                                    values...);
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
               const ParameterVector& values)
      : m_vValues(kNumberOfParameters) {
    if (cov) {
      m_optCovariance = std::move(*cov);
    }
    detail::initialize_parset<parameter_indices_t, params...>::init(*this,
                                                                    values);
  }

  /**
   * @brief copy constructor
   *
   * @param copy object whose content is copied into the new @c ParameterSet
   * object
   */
  ParameterSet(const Self& copy)
      : m_vValues(copy.m_vValues), m_optCovariance(copy.m_optCovariance) {}

  /**
   * @brief move constructor
   *
   * @param copy object whose content is moved into the new @c ParameterSet
   * object
   */
  ParameterSet(Self&& copy) : m_vValues(std::move(copy.m_vValues)) {
    if (copy.m_optCovariance) {
      m_optCovariance = std::move(*copy.m_optCovariance);
    }
  }

  /**
   * @brief standard destructor
   */
  ~ParameterSet() = default;

  /**
   * @brief assignment operator
   *
   * @param rhs object whose content is assigned to this @c ParameterSet object
   */
  Self& operator=(const Self& rhs) {
    m_vValues = rhs.m_vValues;
    m_optCovariance = rhs.m_optCovariance;
    return *this;
  }

  /**
   * @brief move assignment operator
   *
   * @param rhs object whose content is moved into this @c ParameterSet object
   */
  Self& operator=(Self&& rhs) {
    m_vValues = std::move(rhs.m_vValues);
    m_optCovariance = std::move(rhs.m_optCovariance);
    return *this;
  }

  /**
   * @brief swap two objects
   */
  friend void swap(Self& first, Self& second) noexcept {
    using std::swap;
    swap(first.m_vValues, second.m_vValues);
    swap(first.m_optCovariance, second.m_optCovariance);
  }

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
  static constexpr parameter_indices_t getParID() {
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
  ParValue_t getParameter() const {
    return m_vValues(getIndex<parameter>());
  }

  /**
   * @brief access vector with stored parameters
   *
   * @return column vector with @c #kNumberOfParameters rows
   */
  ParameterVector getParameters() const { return m_vValues; }

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
  void setParameter(ParValue_t value) {
    m_vValues(getIndex<parameter>()) =
        ParameterTypeFor<parameter_indices_t, parameter>::type::getValue(value);
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
  void setParameters(const ParameterVector& values) {
    detail::initialize_parset<parameter_indices_t, params...>::init(*this,
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
  ParValue_t getUncertainty() const {
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
   * @brief equality operator
   *
   * @return @c true if stored parameter values are equal and both covariance
   * matrices are
   *         either identical or not set, otherwise @c false
   */
  bool operator==(const Self& rhs) const {
    // shortcut comparison with myself
    if (&rhs == this) {
      return true;
    }

    // parameter values
    if (m_vValues != rhs.m_vValues) {
      return false;
    }
    // both have covariance matrices set
    if ((m_optCovariance && rhs.m_optCovariance) &&
        (*m_optCovariance != *rhs.m_optCovariance)) {
      return false;
    }
    // only one has a covariance matrix set
    if ((m_optCovariance && !rhs.m_optCovariance) ||
        (!m_optCovariance && rhs.m_optCovariance)) {
      return false;
    }

    return true;
  }

  /**
   * @brief inequality operator
   *
   * @return @c true if both objects are not equal, otherwise @c false
   *
   * @sa ParameterSet::operator==
   */
  bool operator!=(const Self& rhs) const { return !(*this == rhs); }

  /**
   * @brief project vector of full parameter set onto parameter sub-space
   *
   * Let \f$ \left(p_1 \dots p_N \right)\f$ be the full set of parameters out of
   * which the \f$m\f$
   * parameters \f$ \left( p_{i_1} \dots p_{i_m} \right), i_1 < i_2 < \dots <
   * i_m, m \le N, i_j \le N\f$
   * are stored in this ParameterSet object. Let \f$ \left(v^0_1 \dots v^0_N
   * \right)\f$ be the parameter
   * values given in the full ParameterSet, then this methods applies the
   * following mapping:
   * \f[
   * \mathbb{R}^{N \times 1} \mapsto \mathbb{R}^{m \times 1} : \left(
   * \begin{array}{c} v_1^0 \\ \vdots \\ v_N^0 \end{array} \right) \mapsto
   * \left( \begin{array}{c} v_{i_1}^0 \\ \vdots \\ v_{i_m}^0 \end{array}
   * \right)
   * \f]
   *
   * @param fullParSet ParameterSet object containing values for all defined
   * parameters
   *
   * @return vector containing only the parameter values from the full parameter
   * vector
   *         which are also defined for this ParameterSet object
   */
  ParameterVector project(const FullParameterSet& fullParSet) const {
    return projector() * fullParSet.getParameters();
  }

  /**
   * @brief calculate residual difference to full parameter vector
   *
   * Calculate the residual differences of the stored parameter values with
   * respect to the corresponding
   * parameter values in the full parameter vector. Hereby, the residual vector
   * is defined as
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
   * vector of parameter values of
   * this ParameterSet object and \f$\vec{v}^0\f$ is the full parameter value
   * vector.
   *
   * @note Constraint and cyclic parameter value ranges are taken into account
   * when calculating
   *       the residual values.
   *
   * @param fullParSet ParameterSet object containing the full set of parameters
   *
   * @return vector containing the residual parameter values of this
   * ParameterSet object
   *         with respect to the given full parameter vector
   *
   * @sa ParameterSet::projector
   */
  /// @cond
  template <
      typename T = Self,
      std::enable_if_t<not std::is_same<T, FullParameterSet>::value, int> = 0>
  /// @endcond
  ParameterVector residual(const FullParameterSet& fullParSet) const {
    return detail::residual_calculator<parameter_indices_t, params...>::result(
        m_vValues, projector() * fullParSet.getParameters());
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
  ParameterVector residual(const Self& otherParSet) const {
    return detail::residual_calculator<parameter_indices_t, params...>::result(
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
  static const ActsMatrix<ParValue_t, kNumberOfParameters, kSizeMax>
  projector() {
    return sProjector;
  }

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
  static void correctValues(ParameterVector& values) {
    detail::value_corrector<params...>::result(values);
  }

 private:
  ParameterVector m_vValues{
      ParameterVector::Zero()};  ///< column vector containing
                                 ///< values of local parameters
  std::optional<CovarianceMatrix> m_optCovariance{
      std::nullopt};  ///< an optional covariance matrix

  static const Projection sProjector;  ///< matrix to project full parameter
                                       /// vector onto local parameter space
};

// initialize static class members
template <typename parameter_indices_t, parameter_indices_t... params>
constexpr unsigned int
    ParameterSet<parameter_indices_t, params...>::kNumberOfParameters;

template <typename parameter_indices_t, parameter_indices_t... params>
const typename ParameterSet<parameter_indices_t, params...>::Projection
    ParameterSet<parameter_indices_t, params...>::sProjector =
        detail::make_projection_matrix<
            kSizeMax, static_cast<unsigned int>(params)...>::init();
}  // namespace Acts
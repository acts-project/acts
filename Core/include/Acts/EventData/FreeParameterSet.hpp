// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// STL include(s)
#include <memory>
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
using FullParameterSet = typename detail::full_freeparset::type;
/// @endcond

/**
 * @class FreeParameterSet
 *
 * @brief Description of a set of (global) parameters
 *
 * @pre
 * The template parameter @c ParameterPolicy must fulfill the following
 * requirements:
 *  -# It must contain a <tt>typedef #FreeID_t</tt> specifying an integral
 * type used to identify different
 *     parameters. This could for example be an @c enum, @c short, or
 * <tt>unsigned int</tt>.
 *     This @c typedef must be convertible to an <tt>unsigned int</tt>
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
 * <tt>Acts::BoundParsDim</tt>.
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
template <FreeID_t... params>
class FreeParameterSet {
 private:
  // local typedefs and constants
  using ParSet_t = FreeParameterSet<params...>;  ///< type of this parameter set
  static constexpr unsigned int NPars =
      sizeof...(params);  ///< number of parameters stored in this class

  // static assert to check that the template parameters are consistent
  static_assert(detail::are_sorted<true, true, FreeID_t, params...>::value,
                "parameter identifiers are not sorted");
  static_assert(
      detail::are_within<unsigned int, 0, BoundParsDim,
                         static_cast<unsigned int>(params)...>::value,
      "parameter identifiers must be greater or "
      "equal to zero and smaller than the total number of parameters");
  static_assert(NPars > 0, "number of stored parameters can not be zero");
  static_assert(
      NPars <= BoundParsDim,
      "number of stored parameters can not exceed number of total parameters");

 public:
  // public typedefs
  /// matrix type for projecting full parameter vector onto local parameter
  /// space
  using Projection_t = ActsMatrix<ParValue_t, NPars, BoundParsDim>;
  /// vector type for stored parameters
  using ParVector_t = ActsVector<ParValue_t, NPars>;
  /// type of covariance matrix
  using CovMatrix_t = ActsSymMatrix<ParValue_t, NPars>;
  /// type for unique pointer to covariance matrix
  using CovPtr_t = std::unique_ptr<const CovMatrix_t>;

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
  FreeParameterSet(CovPtr_t cov,
               std::enable_if_t<sizeof...(Tail) + 1 == NPars, ParValue_t> head,
               Tail... values)
      : m_vValues(NPars), m_pCovariance(std::move(cov)) {
    detail::initialize_parset<FreeID_t, params...>::init(*this, head, values...);
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
  FreeParameterSet(CovPtr_t cov, const ParVector_t& values)
      : m_vValues(NPars), m_pCovariance(std::move(cov)) {
    detail::initialize_parset<FreeID_t, params...>::init(*this, values);
  }

  /**
   * @brief copy constructor
   *
   * @param copy object whose content is copied into the new @c FreeParameterSet
   * object
   */
  FreeParameterSet(const ParSet_t& copy)
      : m_vValues(copy.m_vValues), m_pCovariance(nullptr) {
    if (copy.m_pCovariance) {
      m_pCovariance = std::make_unique<const CovMatrix_t>(*copy.m_pCovariance);
    }
  }

  /**
   * @brief move constructor
   *
   * @param copy object whose content is moved into the new @c FreeParameterSet
   * object
   */
  FreeParameterSet(ParSet_t&& copy)
      : m_vValues(std::move(copy.m_vValues)),
        m_pCovariance(std::move(copy.m_pCovariance)) {}

  /**
   * @brief standard destructor
   */
  ~FreeParameterSet() = default;

  /**
   * @brief assignment operator
   *
   * @param rhs object whose content is assigned to this @c FreeParameterSet object
   */
  ParSet_t& operator=(const ParSet_t& rhs) {
    m_vValues = rhs.m_vValues;
    m_pCovariance =
        (rhs.m_pCovariance
             ? std::make_unique<const CovMatrix_t>(*rhs.m_pCovariance)
             : nullptr);
    return *this;
  }

  /**
   * @brief move assignment operator
   *
   * @param rhs object whose content is moved into this @c FreeParameterSet object
   */
  ParSet_t& operator=(ParSet_t&& rhs) {
    m_vValues = std::move(rhs.m_vValues);
    m_pCovariance = std::move(rhs.m_pCovariance);
    return *this;
  }

  /**
   * @brief swap two objects
   */
  friend void swap(ParSet_t& first, ParSet_t& second) noexcept {
    using std::swap;
    swap(first.m_vValues, second.m_vValues);
    swap(first.m_pCovariance, second.m_pCovariance);
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
  template <FreeID_t parameter>
  static constexpr size_t getIndex() {
    return detail::get_position<FreeID_t, parameter, params...>::value;
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
  static constexpr FreeID_t getParID() {
    return detail::at_index<FreeID_t, index, params...>::value;
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
  template <FreeID_t parameter>
  ParValue_t getParameter() const {
    return m_vValues(getIndex<parameter>());
  }

  /**
   * @brief access vector with stored parameters
   *
   * @return column vector with @c #NPars rows
   */
  ParVector_t getParameters() const { return m_vValues; }

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
  template <FreeID_t parameter>
  void setParameter(ParValue_t value) {
    using parameter_type = typename par_type<parameter>::type;
    m_vValues(getIndex<parameter>()) = parameter_type::getValue(value);
  }

  /**
   * @brief sets values of stored parameters
   *
   * The values of the given vector are interpreted as parameter values in the
   * order
   * of the class template `params...`.
   *
   * @param values vector of length #NPars
   */
  void setParameters(const ParVector_t& values) {
    detail::initialize_parset<FreeID_t, params...>::init(*this, values);
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
  template <FreeID_t parameter>
  bool contains() const {
    return detail::is_contained<FreeID_t, parameter, params...>::value;
  }

  /**
   * @brief access covariance matrix for stored parameters
   *
   * @note The ownership of the covariance matrix is @b not transferred with
   * this call.
   *
   * @return raw pointer to covariance matrix (can be a nullptr)
   */
  const CovMatrix_t* getCovariance() const { return m_pCovariance.get(); }

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
  template <FreeID_t parameter>
  ParValue_t getUncertainty() const {
    if (m_pCovariance) {
      size_t index = getIndex<parameter>();
      return sqrt((*m_pCovariance)(index, index));
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
  void setCovariance(CovPtr_t cov) { m_pCovariance = std::move(cov); }

  /**
   * @brief equality operator
   *
   * @return @c true if stored parameter values are equal and both covariance
   * matrices are
   *         either identical or not set, otherwise @c false
   */
  bool operator==(const ParSet_t& rhs) const {
    // shortcut comparison with myself
    if (&rhs == this) {
      return true;
    }

    // parameter values
    if (m_vValues != rhs.m_vValues) {
      return false;
    }
    // both have covariance matrices set
    if ((m_pCovariance && rhs.m_pCovariance) &&
        (*m_pCovariance != *rhs.m_pCovariance)) {
      return false;
    }
    // only one has a covariance matrix set
    if ((m_pCovariance && !rhs.m_pCovariance) ||
        (!m_pCovariance && rhs.m_pCovariance)) {
      return false;
    }

    return true;
  }

  /**
   * @brief inequality operator
   *
   * @return @c true if both objects are not equal, otherwise @c false
   *
   * @sa FreeParameterSet::operator==
   */
  bool operator!=(const ParSet_t& rhs) const { return !(*this == rhs); }

  /**
   * @brief project vector of full parameter set onto parameter sub-space
   *
   * Let \f$ \left(p_1 \dots p_N \right)\f$ be the full set of parameters out of
   * which the \f$m\f$
   * parameters \f$ \left( p_{i_1} \dots p_{i_m} \right), i_1 < i_2 < \dots <
   * i_m, m \le N, i_j \le N\f$
   * are stored in this FreeParameterSet object. Let \f$ \left(v^0_1 \dots v^0_N
   * \right)\f$ be the parameter
   * values given in the full FreeParameterSet, then this methods applies the
   * following mapping:
   * \f[
   * \mathbb{R}^{N \times 1} \mapsto \mathbb{R}^{m \times 1} : \left(
   * \begin{array}{c} v_1^0 \\ \vdots \\ v_N^0 \end{array} \right) \mapsto
   * \left( \begin{array}{c} v_{i_1}^0 \\ \vdots \\ v_{i_m}^0 \end{array}
   * \right)
   * \f]
   *
   * @param fullParSet FreeParameterSet object containing values for all defined
   * parameters
   *
   * @return vector containing only the parameter values from the full parameter
   * vector
   *         which are also defined for this FreeParameterSet object
   */
  ParVector_t project(const FullFreeParameterSet& fullParSet) const {
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
   * this FreeParameterSet object and \f$\vec{v}^0\f$ is the full parameter value
   * vector.
   *
   * @note Constraint and cyclic parameter value ranges are taken into account
   * when calculating
   *       the residual values.
   *
   * @param fullParSet FreeParameterSet object containing the full set of parameters
   *
   * @return vector containing the residual parameter values of this
   * FreeParameterSet object
   *         with respect to the given full parameter vector
   *
   * @sa FreeParameterSet::projector
   */
  /// @cond
  template <
      typename T = ParSet_t,
      std::enable_if_t<not std::is_same<T, FullFreeParameterSet>::value, int> = 0>
  /// @endcond
  ParVector_t residual(const FullFreeParameterSet& fullParSet) const {
    return detail::residual_calculator<params...>::result(
        m_vValues, projector() * fullParSet.getParameters());
  }

  /**
   * @brief calculate residual difference to other parameter vector
   *
   * Calculate the residual differences of the stored parameter values with
   * respect to the values of
   * another FreeParameterSet object containing the same set of parameters. Hereby,
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
   * where \f$\vec{v}\f$ is the vector of parameter values of this FreeParameterSet
   * object and \f$\vec{v}^0\f$
   * is the parameter value vector of the other FreeParameterSet object.
   *
   * @note Constraint and cyclic parameter value ranges are taken into account
   * when calculating
   *       the residual values.
   *
   * @param otherParSet FreeParameterSet object with identical set of contained
   * parameters
   *
   * @return vector containing the residual parameter values of this
   * FreeParameterSet object
   *         with respect to the given other parameter set
   */
  ParVector_t residual(const ParSet_t& otherParSet) const {
    return detail::residual_calculator<params...>::result(
        m_vValues, otherParSet.m_vValues);
  }

  /**
   * @brief get projection matrix
   *
   * The projection matrix performs a mapping of the full parameter space onto
   * the sub-space
   * spanned by the parameters defined in this FreeParameterSet object.
   *
   * @return constant matrix with @c #NPars rows and @c #Acts::BoundParsDim
   * columns
   */
  static const ActsMatrix<ParValue_t, NPars, BoundParsDim> projector() {
    return sProjector;
  }

  /**
   * @brief number of stored parameters
   *
   * @return number of stored parameters
   */
  static constexpr unsigned int size() { return NPars; }

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
  static void correctValues(ParVector_t& values) {
    detail::value_corrector<params...>::result(values);
  }

 private:
  ParVector_t
      m_vValues;  ///< column vector containing values of local parameters
  CovPtr_t m_pCovariance;  ///< unique pointer to covariance matrix

  static const Projection_t sProjector;  ///< matrix to project full parameter
                                         /// vector onto local parameter space
};

// initialize static class members
template <FreeID_t... params>
constexpr unsigned int FreeParameterSet<params...>::NPars;

template <FreeID_t... params>
const typename FreeParameterSet<params...>::Projection_t
    FreeParameterSet<params...>::sProjector = detail::make_projection_matrix<
        BoundParsDim, static_cast<unsigned int>(params)...>::init();
}  // namespace Acts
#ifndef ACTS_PARAMETERSET_H
#define ACTS_PARAMETERSET_H 1

// STL include(s)
#include <memory>
#include <type_traits>
#include <utility>

#include "ACTS/Utilities/Definitions.hpp"
// ACTS includes
#include "ACTS/Utilities/ParameterDefinitions.hpp"

namespace Acts
{
  /// @cond
  // forward declaration
  template<ParID_t... params>
  class ParameterSet;
  /// @endcond

  /// @cond DEV
  /**
   * @brief anonymous namespace for implementation details
   *
   * Anonymous namespace containing implementation details for consistency checks at compile time.
   */
  namespace
  {
    /**
     * @brief check whether integral values are sorted
     *
     * @tparam ascending boolean flag to check for ascending order (@c true) or descending order (@c false)
     * @tparam strict boolean flag whether strict ordering is required
     * @tparam T integral type of values whose order should be checked
     * @tparam values template parameter pack containing the list of values
     *
     * @test Unit tests are implemented \link Acts::Test::BOOST_AUTO_TEST_CASE(are_sorted_helper_tests) here\endlink.
     *
     * @return `are_sorted<asc,strict,T,values...>::value` is @c true if the given values are properly sorted,
     *         otherwise @c false
     */
    template<bool ascending, bool strict, typename T, T ... values>
    struct are_sorted;

    /**
     * @brief check whether integral values are within a given range
     *
     * @tparam T integral type of values whose range should be checked
     * @tparam MIN lower accepted bound of values (inclusive)
     * @tparam MAX upper accepted bound of values (exclusive)
     * @tparam values template parameter pack containing the list of values
     *
     * @test Unit tests are implemented \link Acts::Test::BOOST_AUTO_TEST_CASE(are_within_helper_tests) here\endlink.
     *
     * @return `are_within<T,MIN,MAX,values...>::value` is @c true if all given values are within the
     *          interval [MIN,MAX), otherwise @c false
     */
    template<typename T, T MIN, T MAX, T... values>
    struct are_within;

    /**
     * @brief generate ParameterSet type containing all defined parameters
     *
     * @return `full_parset<Policy>::type` is equivalent to `ParameterSet<Policy,ID_t(0),ID_t(1),...,ID_t(N-1)>`
     *         where @c ID_t is a @c typedef to `Policy::par_id_type` and @c N is the total number of parameters
     */
    struct full_parset
    {
      template<ParID_t v,typename C>
      struct add_to_value_container;

      template<ParID_t v,ParID_t... others>
      struct add_to_value_container<v,std::integer_sequence<ParID_t,others...> >
      {
        typedef std::integer_sequence<ParID_t,others...,v> type;
      };

      template<typename T,unsigned int N>
      struct tparam_generator
      {
        typedef typename add_to_value_container<static_cast<T>(N),typename tparam_generator<T,N-1>::type>::type type;
      };

      template<typename T>
      struct tparam_generator<T,0>
      {
        typedef std::integer_sequence<T,static_cast<T>(0)> type;
      };

      template<typename T>
      struct converter;

      template<ParID_t... values>
      struct converter<std::integer_sequence<ParID_t,values...> >
      {
        typedef ParameterSet<values...> type;
      };

      typedef typename converter<typename tparam_generator<ParID_t,Acts::NGlobalPars-1>::type>::type type;
    };
    /// @endcond
  }  // end of anonymous namespace
  /// @endcond

  typedef typename full_parset::type FullParameterSet;

  /**
   * @class ParameterSet
   *
   * @brief Description of a set of (local) parameters
   *
   * @pre
   * The template parameter @c ParameterPolicy must fulfill the following requirements:
   *  -# It must contain a <tt>typedef #par_id_type</tt> specifying an integral type used to identify different
   *     parameters. This could for example be an @c enum, @c short, or <tt>unsigned int</tt>.
   *     This @c typedef must be convertible to an <tt>unsigned int</tt>
   *  -# It must contain a <tt>typedef #par_value_type</tt> specifying the type of the parameter values. This could for
   *     instance be @c double, or @c float.
   *  -# It must contain a definition of an integral constant named @c N which is assignable to an <tt>unsigned
   *     int</tt> and which is equal to the total number of parameters in the system.
   * @pre
   *
   * The template parameter pack @c params must be given in a strictly ascending order. The parameter pack must
   * be non-empty and it cannot contain more elements than <tt>Acts::NGlobalPars</tt>.
   *
   * @test The behavior of this class is tested in the following unit tests:
   *       - \link Acts::Test::BOOST_AUTO_TEST_CASE(parset_consistency_tests) general consistency\endlink
   *       - \link Acts::Test::BOOST_AUTO_TEST_CASE(parset_copy_assignment_tests) copy/assignment/swap\endlink
   *       - \link Acts::Test::BOOST_AUTO_TEST_CASE(parset_comparison_tests) comparison operators\endlink
   *       - \link Acts::Test::BOOST_AUTO_TEST_CASE(parset_projection_tests) projection matrices\endlink
   *       - \link Acts::Test::BOOST_AUTO_TEST_CASE(parset_residual_tests) residual calculation\endlink
   *
   * @tparam ParameterPolicy  struct or class containing the parameter definitions (see above)
   * @tparam params           parameter pack containing the (local) parameters stored in this class
   */
  template<ParID_t ... params>
  class ParameterSet
  {
  private:
    // local typedefs and constants
    typedef ParameterSet<params...> ParSet_t;        ///< type of this parameter set
    static constexpr unsigned int NPars = sizeof...(params);         ///< number of parameters stored in this class

    // static assert to check that the template parameter are consistent
    static_assert(are_sorted<true,true,ParID_t,params...>::value,"parameter identifiers are not sorted");
    static_assert(are_within<unsigned int,0,Acts::NGlobalPars,static_cast<unsigned int>(params)...>::value,"parameter identifiers must be greater or "
        "equal to zero and smaller than the total number of parameters");
    static_assert(NPars > 0,"number of stored parameters can not be zero");
    static_assert(NPars <= Acts::NGlobalPars,"number of stored parameters can not exceed number of total parameters");

  public:
    // public typedefs
    typedef ActsMatrix<ParValue_t,NPars,Acts::NGlobalPars> Projection_t;       ///< matrix type for projecting full parameter vector onto local parameter space
    typedef ActsVector<ParValue_t,NPars> ParVector_t;                 ///< vector type for stored parameters
    typedef ActsSymMatrix<ParValue_t,NPars> CovMatrix_t;              ///< type of covariance matrix
    typedef std::unique_ptr<CovMatrix_t> Cov_uptr;                   ///< type for unique pointer to covariance matrix

    /**
     * @brief initialize values of stored parameters and their covariance matrix
     *
     * @note  No validation of the given covariance matrix is performed.
     *
     * @param cov unique pointer to covariance matrix (nullptr is accepted)
     * @param values parameter pack with values for the stored parameters
     */
    template<typename ... Tail>
    ParameterSet(Cov_uptr cov,std::enable_if_t<sizeof...(Tail) + 1 == NPars, ParValue_t> head, Tail ... values);

    /**
     * @brief initialize parameter values from vector and set their covariance matrix
     *
     * @note The values in the passed vector are interpreted as parameter values in the order given
     *       by the class template @c params. No validation of the given covariance matrix is performed.
     *
     * @param cov unique pointer to covariance matrix (nullptr is accepted)
     * @param values vector with parameter values
     */
    ParameterSet(Cov_uptr cov,const ParVector_t& values);

    /**
     * @brief copy constructor
     *
     * @param copy object whose content is copied into the new @c ParameterSet object
     */
    ParameterSet(const ParSet_t& copy);

    /**
     * @brief move constructor
     *
     * @param copy object whose content is moved into the new @c ParameterSet object
     */
    ParameterSet(ParSet_t&& copy);

    /**
     * @brief standard destructor
     */
    ~ParameterSet() = default;

    /**
     * @brief assignment operator
     *
     * @param rhs object whose content is assigned to this @c ParameterSet object
     */
    ParSet_t& operator=(const ParSet_t& rhs);

    /**
     * @brief move assignment operator
     *
     * @param rhs object whose content is moved into this @c ParameterSet object
     */
    ParSet_t& operator=(ParSet_t&& rhs);

    /**
     * @brief swap two objects
     */
    friend void swap(ParSet_t& first,ParSet_t& second) noexcept
    {
      using std::swap;
      swap(first.m_vValues,second.m_vValues);
      swap(first.m_pCovariance,second.m_pCovariance);
    }

    /**
     * @brief retrieve stored value for given parameter
     *
     * @tparam parameter identifier for the parameter to be retrieved
     * @remark @c parameter must be part of the template parameter pack @c params. Otherwise a compile-time
     *         error is generated.
     *
     * @return value of the stored parameter
     */
    template<ParID_t parameter>
    ParValue_t getParameter() const;

    /**
     * @brief access vector with stored parameters
     *
     * @return column vector with @c #NPars rows
     */
    ParVector_t getParameters() const;

    /**
     * @brief sets value for given parameter
     *
     * @tparam parameter identifier for the parameter to be stored
     * @remark @c parameter must be part of the template parameter pack @c params. Otherwise a compile-time
     *         error is generated.
     *
     * @return previously stored value of this parameter
     */
    template<ParID_t parameter>
    void setParameter(ParValue_t value);

    /**
     * @brief sets values of stored parameters
     *
     * The values of the given vector are interpreted as parameter values in the order
     * of the class template `params...`.
     *
     * @param values vector of length #NPars
     */
    void setParameters(const ParVector_t& values);

    /**
     * @brief checks whether a given parameter is included in this set of parameters

     * @tparam parameter identifier for the parameter to be retrieved
     * @remark @c parameter must be part of the template parameter pack @c params. Otherwise a compile-time
     *         error is generated.
     *
     * @return @c true if the parameter is stored in this set, otherwise @c false
     */
    template<ParID_t parameter>
    bool contains() const;

    /**
     * @brief access covariance matrix for stored parameters
     *
     * @note The ownership of the covariance matrix is @b not transferred with this call.
     *
     * @return raw pointer to covariance matrix (can be zero)
     */
    const CovMatrix_t* getCovariance() const;

    /**
     * @brief access uncertainty for individual parameter
     *
     * @tparam parameter identifier for the parameter to be retrieved
     * @remark @c parameter must be part of the template parameter pack @c params. Otherwise a compile-time
     *         error is generated.
     *
     * @return uncertainty \f$\sigma \ge 0\f$ of given parameter, a negative value is returned if no
     *         covariance matrix is set
     */
    template<ParID_t parameter>
    ParValue_t uncertainty() const;

    /**
     * @brief update covariance matrix
     *
     * @note No validation of the given covariance matrix is performed.
     *
     * @param cov unique pointer to new covariance matrix (nullptr is accepted)
     */
    void setCovariance(Cov_uptr cov);

    /**
     * @brief equality operator
     *
     * @return @c true if stored parameter values are equal and both covariance matrices are
     *         either identical or not set, otherwise @c false
     */
    bool operator==(const ParSet_t& rhs) const;

    /**
     * @brief inequality operator
     *
     * @return @c true if both objects are not equal, otherwise @c false
     *
     * @sa ParameterSet::operator==
     */
    bool operator!=(const ParSet_t& rhs) const
    {
      return !(*this == rhs);
    }

    /**
     * @brief project vector of full parameter set onto parameter sub-space
     *
     * Let \f$ \left(p_1 \dots p_N \right)\f$ be the full set of parameters out of which the \f$m\f$
     * parameters \f$ \left( p_{i_1} \dots p_{i_m} \right), i_1 < i_2 < \dots < i_m, m \le N, i_j \le N\f$
     * are stored in this ParameterSet object. Let \f$ \left(v^0_1 \dots v^0_N \right)\f$ be the parameter
     * values given in the full ParameterSet, then this methods applies the following mapping:
     * \f[
     * \mathbb{R}^{N \times 1} \mapsto \mathbb{R}^{m \times 1} : \left( \begin{array}{c} v_1^0 \\ \vdots \\ v_N^0 \end{array} \right) \mapsto \left( \begin{array}{c} v_{i_1}^0 \\ \vdots \\ v_{i_m}^0 \end{array} \right)
     * \f]
     *
     * @param fullParSet ParameterSet object containing values for all defined parameters
     *
     * @return vector containing only the parameter values from the full parameter vector
     *         which are also defined for this ParameterSet object
     */
    ParVector_t project(const FullParameterSet& fullParSet) const;

    /**
     * @brief calculate residual difference to full parameter vector
     *
     * Calculate the residual differences of the stored parameter values with respect to the corresponding
     * parameter values in the full parameter vector. Hereby, the residual vector is defined as
     *
     * \f[
     * \vec{r} = \left( \begin{array}{c} r_{i_1} \\ \vdots \\ r_{i_m} \end{array} \right)
     *  = \left( \begin{array}{c} v_{i_1} \\ \vdots \\ v_{i_m} \end{array} \right) -  \mathrm{Proj} \left( \begin{array}{c} v^0_{1} \\ \vdots \\ v^0_{N} \end{array} \right)
     *  = \vec{v} - \mathrm{Proj} \left( \vec{v}^0 \right)
     * \f]
     *
     * where \f$\mathrm{Proj}\f$ is the projection matrix, \f$\vec{v}\f$ is the vector of parameter values of
     * this ParameterSet object and \f$\vec{v}^0\f$ is the full parameter value vector.
     *
     * @note Constraint and cyclic parameter value ranges are taken into account when calculating
     *       the residual values.
     *
     * @param fullParSet ParameterSet object containing the full set of parameters
     *
     * @return vector containing the residual parameter values of this ParameterSet object
     *         with respect to the given full parameter vector
     *
     * @sa ParameterSet::projector
     */
    /// @cond
    template<typename T = ParSet_t,std::enable_if_t<not std::is_same<T,FullParameterSet>::value,int> = 0>
    /// @endcond
    ParVector_t residual(const FullParameterSet& fullParSet) const;

    /**
     * @brief calculate residual difference to other parameter vector
     *
     * Calculate the residual differences of the stored parameter values with respect to the values of
     * another ParameterSet object containing the same set of parameters. Hereby, the residual vector is
     * defined as
     *
     * \f[
     * \vec{r} = \left( \begin{array}{c} r_{i_1} \\ \vdots \\ r_{i_m} \end{array} \right)
     *  = \left( \begin{array}{c} v_{i_1} \\ \vdots \\ v_{i_m} \end{array} \right) -  \left( \begin{array}{c} v^0_{1} \\ \vdots \\ v^0_{N} \end{array} \right)
     *  = \vec{v} - \left( \vec{v}^0 \right)
     * \f]
     *
     * where \f$\vec{v}\f$ is the vector of parameter values of this ParameterSet object and \f$\vec{v}^0\f$
     * is the parameter value vector of the other ParameterSet object.
     *
     * @note Constraint and cyclic parameter value ranges are taken into account when calculating
     *       the residual values.
     *
     * @param otherParSet ParameterSet object with identical set of contained parameters
     *
     * @return vector containing the residual parameter values of this ParameterSet object
     *         with respect to the given other parameter set
     */
    ParVector_t residual(const ParSet_t& otherParSet) const;

    /**
     * @brief get projection matrix
     *
     * The projection matrix performs a mapping of the full parameter space onto the sub-space
     * spanned by the parameters defined in this ParameterSet object.
     *
     * @return constant matrix with @c #NPars rows and @c #Acts::NGlobalPars columns
     */
    static const ActsMatrix<ParValue_t,NPars,Acts::NGlobalPars> projector()
    {
      return sProjector;
    }

    /**
     * @brief number of stored parameters
     *
     * @return number of stored parameters
     */
    static constexpr unsigned int size()
    {
      return NPars;
    }

    /**
     * @brief correct given parameter values
     *
     * Check that the given values are within in a valid range for the corresponding parameter. If not, an
     * in-place correction is applied. The values are interpreted as parameter values in the same order as
     * specified in the class template @c params.
     *
     * @param values vector with parameter values to be checked and corrected if necessary
     */
    static void correctValues(ParVector_t& values);

  private:
    ParVector_t m_vValues;           ///< column vector containing values of local parameters
    Cov_uptr    m_pCovariance;       ///< unique pointer to covariance matrix

    static const Projection_t sProjector;  ///< matrix to project full parameter vector onto local parameter space
  };
} // end of namespace Acts

#include "ACTS/EventData/detail/ParameterSet.icc"

#endif // ACTS_PARAMETERSET_H

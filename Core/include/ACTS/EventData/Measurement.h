#ifndef ACTS_MEASUREMENT_H
#define ACTS_MEASUREMENT_H 1

// STL include(s)
#include <type_traits>
#include <utility>
#include <memory>

// boost include(s)
#include <boost/mpl/vector.hpp>

// ACTS includes
#include "CoreUtils/ParameterDefinitions.h"
#include "ParameterSet/ParameterSet.h"
#include "TrackParameters/TrackParameters.h"

// boost include(s)
#include <boost/variant.hpp>

namespace Acts
{
  // forward declarations
  class Surface;
}

namespace Acts
{
  /**
   * @brief base class for Measurements
   *
   * This class describes the measurement of track parameters at a certain Surface in the
   * TrackingGeometry.
   *
   * @test The behaviour of this class is tested in the following unit test:
   *       - \link Acts::Test::BOOST_AUTO_TEST_CASE(measurement_initialization) initialization\endlink
   * @tparam Identifier identification object for this measurement
   * @tparam params     parameter pack containing the measured parameters
   */
  template<typename Identifier,ParID_t... params>
  class Measurement
  {
  private:
    // private typedef's
    typedef ParameterSet<params...>    ParSet_t;   ///< type of the underlying ParameterSet object

  public:
    typedef typename ParSet_t::ParVector_t ParVector_t;            ///< type of the vector containing the parameter values
    typedef typename ParSet_t::CovMatrix_t CovMatrix_t;            ///< type of the covariance matrix of the measurement

    /**
     * @brief standard constructor
     *
     * Interface class for all possible measurements.
     *
     * @note The given ParameterSet object is copied while only a reference to the given surface is stored.
     *       The user must ensure that the lifetime of the @c Surface object surpasses the lifetime of this Measurement
     *       object.<br />
     *       The covariance matrix object is moved into this Measurement object and should not be used afterwards.<br />
     *       The given parameter values are interpreted as values to the parameters as defined in the class template
     *       argument @c params.
     *
     * @attention The current design will fail if the in-memory location of the @c Surface object is changed (e.g.
     *            if it is stored in a container and this gets relocated).
     *
     * @param surface surface at which the measurement took place
     * @param cov covariance matrix of the measurement.
     * @param head,values consistent number of parameter values of the measurement
     */
    template<typename ... Tail>
    Measurement(const Surface& surface,
                const Identifier& id,
                CovMatrix_t&& cov,
                typename std::enable_if<sizeof...(Tail) + 1 == sizeof...(params),ParValue_t>::type head,
                Tail... values):
      m_oParameters(std::make_unique<CovMatrix_t>(std::move(cov)),head,values...),
      m_pSurface(&surface),
      m_oIdentifier(id)
    {}


    /**
     * @brief virtual destructor
     */
    virtual ~Measurement() = default;

    /**
     * @brief copy constructor
     */
    Measurement(const Measurement<Identifier,params...>& copy):
      m_oParameters(copy.m_oParameters),
      m_pSurface(copy.m_pSurface),
      m_oIdentifier(copy.m_oIdentifier)
    {}

    /**
     * @brief move constructor
     */
    Measurement(Measurement<Identifier,params...>&& rhs):
      m_oParameters(std::move(rhs.m_oParameters)),
      m_pSurface(rhs.m_pSurface),
      m_oIdentifier(std::move(rhs.m_oIdentifier))
    {}

    /**
     * @brief copy assignment operator
     */
    Measurement<Identifier,params...>& operator=(const Measurement<Identifier,params...>& rhs)
    {
      m_oParameters = rhs.m_oParameters;
      m_pSurface    = rhs.m_pSurface;
      m_oIdentifier = rhs.m_oIdentifier;

      return *this;
    }

    /**
     * @brief move assignment operator
     */
    Measurement<Identifier,params...>& operator=(Measurement<Identifier,params...>&& rhs)
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
     * @remark @c parameter must be part of the template parameter pack @c params. Otherwise a compile-time
     *         error is generated.
     *
     * @return value of the stored parameter
     */
    template<ParID_t parameter>
    ParValue_t get() const
    {
      return m_oParameters.template getParameter<parameter>();
    }

    /**
     * @brief access vector with measured parameter values
     *
     * @return column vector whose size is equal to the dimensionality of this Measurement. The values are
     *         given for the measured parameters in the order defined by the class template argument @c params.
     */
    ParVector_t parameters() const
    {
      return m_oParameters.getParameters();
    }

    /**
     * @brief access covariance matrix of the measured parameter values
     *
     * @pre The stored ParameterSet object must have a valid pointer to covariance matrix assigned.
     *
     * @return covariance matrix of the measurement
     */
    CovMatrix_t covariance() const
    {
      return *m_oParameters.getCovariance();
    }

    /**
     * @brief number of measured parameters
     *
     * @return number of measured parameters
     */
    static constexpr unsigned int size()
    {
      return ParSet_t::size();
    }

    /**
     * @brief access associated surface
     *
     * @pre The @c Surface object used to construct this @c Measurement object must still be valid
     *      at the same memory location.
     *
     * @return reference to surface at which the measurement took place
     */
    const Acts::Surface& associatedSurface() const
    {
      return *m_pSurface;
    }

    /**
     * @brief access to global measurement identifier
     *
     * @return identifier object
     */
    Identifier identifier() const
    {
      return m_oIdentifier;
    }

    /**
     * @brief calculate residual with respect to given track parameters
     *
     * @note It is checked that the residual for non-local parameters are in valid range (e.g.
     *       residuals in \f$\phi\f$ are corrected).
     *
     * @todo Implement check that TrackParameters are defined at the same Surface as the Measurement is.
     * @todo Implement validity check for residuals of local parameters.
     *
     * @param trackPars reference TrackParameters object
     *
     * @return vector with the residual parameter values (in valid range)
     *
     * @sa ParameterSet::residual
     */
    ParVector_t residual(const TrackParameters& trackPars) const
    {
      return m_oParameters.residual(trackPars.getParameterSet());
    }

    /**
     * @brief equality operator
     *
     * @return @c true if parameter sets and associated surfaces compare equal, otherwise @c false
     */
    virtual bool operator==(const Measurement<Identifier,params...>& rhs) const
    {
      return ((m_oParameters == rhs.m_oParameters) &&
              (*m_pSurface == *rhs.m_pSurface) &&
              (m_oIdentifier == rhs.m_oIdentifier));
    }

    /**
     * @brief inequality operator
     *
     * @return @c true if both objects are not equal, otherwise @c false
     *
     * @sa Measurement::operator==
     */
    bool operator!=(const Measurement<Identifier,params...>& rhs) const
    {
      return !(*this == rhs);
    }

  private:
    ParSet_t       m_oParameters;   ///< measured parameter set
    const Surface* m_pSurface;      ///< surface at which the measurement took place
    Identifier     m_oIdentifier;   ///< identifier for this measurement
  };

  /// @cond DEV
  namespace detail
  {
    /**
     * @brief generate boost::variant type for all possible Measurement's
     */
    template<typename ID>
    struct fittable_type_generator;

    /// @cond
    template<typename ID>
    struct fittable_type_generator
    {
      template<ParID_t... params>
      using Meas_t = Measurement<ID,params...>;

      template<typename... T>
      struct container
      {};

      template<typename T,typename U>
      struct add_prepended;

      template<ParID_t first,typename... others>
      struct add_prepended<Meas_t<first>,container<others...> >
      {
        typedef container<typename add_prepended<Meas_t<first>,others>::type...,others...> type;
      };

      template<ParID_t first,ParID_t... others>
      struct add_prepended<Meas_t<first>,Meas_t<others...> >
      {
        typedef Meas_t<first,others...> type;
      };

      template<ParID_t... first>
      struct add_prepended<Meas_t<first...>,boost::mpl::na>
      {
        typedef Meas_t<first...> type;
      };

      template<typename T,typename C>
      struct add_to_container;

      template<typename T,typename... others>
      struct add_to_container<T,container<others...> >
      {
        typedef container<T,others...> type;
      };

      template<typename T>
      struct generator_impl;

      template<ParID_t first,ParID_t... others>
      struct generator_impl<container<Meas_t<first>, Meas_t<others...> > >
      {
        typedef container<Meas_t<first>,Meas_t<others...>,Meas_t<first,others...> > type;
      };

      template<ParID_t first,typename next,typename... others>
      struct generator_impl<container<Meas_t<first>, next, others...> >
      {
        typedef typename generator_impl<container<next, others...> >::type others_combined;
        typedef typename add_prepended<Meas_t<first>,others_combined>::type prepended;
        typedef typename add_to_container<Meas_t<first>,prepended>::type type;
      };

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
        typedef typename add_to_value_container<static_cast<ParID_t>(N),typename tparam_generator<T,N-1>::type>::type type;
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
        typedef container<Meas_t<values>...> type;
      };

      template<typename... types>
      struct to_boost_vector;

      template<typename first,typename... rest>
      struct to_boost_vector<first,rest...>
      {
        typedef typename boost::mpl::push_front<typename to_boost_vector<rest...>::type,first>::type type;
      };

      template<typename last>
      struct to_boost_vector<last>
      {
        typedef boost::mpl::vector<last> type;
      };

      template<typename... MeasTypes>
      struct converter<container<MeasTypes...> >
      {

        typedef typename boost::make_variant_over<typename to_boost_vector<MeasTypes...>::type>::type type;
      };

      typedef typename tparam_generator<ParID_t,Acts::NGlobalPars-1>::type par_list;
      typedef typename converter<par_list>::type meas_list;
      typedef typename generator_impl<meas_list>::type permutations;
      typedef typename converter<permutations>::type type;
    };
    /// @endcond
  } // end of namespace details
  /// @endcond

  /**
   * @brief general type for any possible Measurement
   */
  template<typename Identifier>
  using FittableMeasurement = typename detail::fittable_type_generator<Identifier>::type;
} // end of namespace Acts

#endif // ACTS_MEASUREMENT_H



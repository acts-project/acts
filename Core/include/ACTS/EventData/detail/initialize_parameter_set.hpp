#ifndef ACTS_INITIALIZE_PARAMETER_SET_H
#define ACTS_INITIALIZE_PARAMETER_SET_H 1

namespace Acts
{
  /// @cond detail
  namespace detail
  {
    /**
     * @brief initialize parameter set with given parameter values
     *
     * @note uses the templated ParameterSet::set method for assigning the individual components
     *
     * Possible invocations are:
     * * `initialize<T,params...>::init(parSet,values...)` where `parSet` is the ParameterSet object to be
     *    initialized and `values` are a consistent number of parameter values (with compatible type)
     * * `initialize<T,params...>::init(parSet,values)` where `parSet` is the ParameterSet object to be
     *    initialized and `values` is an Eigen vector of consistent size
     *
     * @tparam T type of the parameters stored in the corresponding @c ParameterSet class
     * @tparam params template parameter pack containing the multiple identifiers
     */
    template<typename T, T... params>
    struct initialize_parset;

    /// @cond
    template<typename T,T first,T... others>
    struct initialize_parset<T,first,others...>
    {
      template<typename ParSetType,typename first_value_type,typename... other_value_types>
      static void init(ParSetType& parSet,const first_value_type& v1, const other_value_types&... values)
      {
        parSet.template setParameter<first>(v1);
        initialize_parset<T,others...>::init(parSet,values...);
      }

      template<typename ParSetType>
      static void init(ParSetType& parSet,const typename ParSetType::ParVector_t& values, const unsigned int& pos = 0)
      {
        parSet.template setParameter<first>(values(pos));
        initialize_parset<T,others...>::init(parSet,values,pos+1);
      }
    };

    template<typename T,T last>
    struct initialize_parset<T,last>
    {
      template<typename ParSet_tType,typename last_value_type>
      static void init(ParSet_tType& ParSet_t,const last_value_type& v1)
      {
        ParSet_t.template setParameter<last>(v1);
      }

      template<typename ParSetType>
      static void init(ParSetType& parSet,const typename ParSetType::ParVector_t& values, const unsigned int& pos = 0)
      {
        parSet.template setParameter<last>(values(pos));
      }
    };
    /// @endcond
  }  // end of namespace detail
  /// @endcond
}  // end of namespace Acts
#endif // ACTS_INITIALIZE_PARAMETER_SET_H

#ifndef ACTS_RESIDUAL_CALCULATOR_H
#define ACTS_RESIDUAL_CALCULATOR_H 1

// ACTS include(s)
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"

namespace Acts
{
  /// @cond detail
  namespace detail
  {
    /**
     * @brief calculate residuals from two parameter vectors
     *
     * Calculate the difference between the two given vectors with parameter values. Possible
     * corrections for bounded or cyclic parameters are applied.
     *
     * @tparam params template parameter pack containing the multiple parameter identifiers
     *
     * @return `residual_calculator<params...>::result(first,second)` yields the residuals of
     *         `first` with respect to `second`
     */
    template<ParID_t... params>
    struct residual_calculator;

    /// @cond
    template<typename R,ParID_t... params>
    struct residual_calculator_impl;

    template<ParID_t... params>
    struct residual_calculator
    {
      typedef ActsVector<ParValue_t,sizeof...(params)> ParVector_t;

      static ParVector_t result(const ParVector_t& test,const ParVector_t& ref)
      {
        ParVector_t result;
        residual_calculator_impl<ParVector_t,params...>::calculate(result,test,ref,0);
        return result;
      }
    };

    template<typename R, ParID_t first,ParID_t... others>
    struct residual_calculator_impl<R,first,others...>
    {
      static void calculate(R& result,const R& test,const R& ref,unsigned int pos)
      {
        typedef typename par_type<first>::type parameter_type;
        result(pos) = parameter_type::getDifference(test(pos),ref(pos));
        residual_calculator_impl<R,others...>::calculate(result,test,ref,pos+1);
      }
    };

    template<typename R,ParID_t last>
    struct residual_calculator_impl<R,last>
    {
      static void calculate(R& result,const R& test,const R& ref,unsigned int pos)
      {
        typedef typename par_type<last>::type parameter_type;
        result(pos) = parameter_type::getDifference(test(pos),ref(pos));
      }
    };
    /// @endcond
  }  // end of namespace detail
  /// @endcond
}  // end of namespace Acts

#endif // ACTS_RESIDUAL_CALCULATOR_H

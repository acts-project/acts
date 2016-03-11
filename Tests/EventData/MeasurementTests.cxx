/**
 * @file MeasurementTests.cxx
 */

// Boost include(s)
#define BOOST_TEST_MODULE Measurement Tests
#include <boost/test/included/unit_test.hpp>

// ATS include(s)
#include "CoreUtils/ParameterDefinitions.h"
#include "Measurement/Measurement.h"
#include "Surfaces/CylinderSurface.h"

namespace Ats
{
  namespace Test
  {
    template<ParID_t... params>
    using Measurement_t = Measurement<unsigned long int,params...>;

    /**
     * @brief Unit test for creation of Measurement object
     */
    BOOST_AUTO_TEST_CASE(measurement_initialization)
    {
      CylinderSurface cylinder(3,10);
      AtsSymMatrixD<2> cov;
      cov << 0.04,0,
             0,0.1;
      Measurement_t<ParDef::eLOC_1,ParDef::eLOC_2> m(cylinder,0,std::move(cov),-0.1,0.45);
    }
  }  // end of namespace Test
}  // end of namespace Ats

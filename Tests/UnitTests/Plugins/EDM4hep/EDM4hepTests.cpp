#include <boost/test/unit_test.hpp>

#include "Acts/Plugins/EDM4hep/test.hpp"

using namespace Acts;

BOOST_AUTO_TEST_SUITE(EDM4hep)

BOOST_AUTO_TEST_CASE(Tests) {
  open("/home/andreas/cern/source/OpenDataDetector/output_edm4hep.root");
}

BOOST_AUTO_TEST_SUITE_END()

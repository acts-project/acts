// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE Material Tests
#include <boost/test/included/unit_test.hpp>
#include <climits>
#include "Acts/Material/Material.hpp"

namespace Acts {

namespace Test {

  // the maximum tolerance is half the accuracy
  float elMaxTolerance = 0.5 / float(UCHAR_MAX);

  BOOST_AUTO_TEST_CASE(ElementFraction_test)
  {
    // carbon parameters, atomic charge is Z
    unsigned char carbonZ = 12;
    // a fraction between 0 and 255
    unsigned char carbonFractionCH    = 46;
    float         carbonFractionFloat = float(46. / UCHAR_MAX);

    // test the carbon fraction
    ElementFraction carbonEFC(carbonZ, carbonFractionCH);
    ElementFraction carbonEFF((unsigned int)carbonZ, carbonFractionFloat);

    // check if you get the element and the fraction back
    BOOST_CHECK_EQUAL(12ul, carbonEFC.element());
    BOOST_CHECK_CLOSE(
        carbonFractionFloat, carbonEFC.fraction(), elMaxTolerance);
    BOOST_CHECK_EQUAL(12ul, carbonEFF.element());
    BOOST_CHECK_CLOSE(
        carbonFractionFloat, carbonEFF.fraction(), elMaxTolerance);
  }

  BOOST_AUTO_TEST_CASE(MaterialComposition_test)
  {

    // Carbon fraction
    unsigned int    carbonZ        = 12;
    float           carbonFraction = 0.45;
    ElementFraction carbon(carbonZ, carbonFraction);

    // Silicon fraction
    unsigned int    siliconZ       = 14;
    float           siliconFracton = 0.1;
    ElementFraction silicon(siliconZ, siliconFracton);

    // Titanium fraction
    unsigned int    titantiumZ       = 22;
    float           titaniumFraction = 0.25;
    ElementFraction titanium(titantiumZ, titaniumFraction);

    // Copper fraction
    unsigned int    copperZ        = 29;
    float           copperFraction = 0.2;
    ElementFraction copper(copperZ, copperFraction);

    std::vector<ElementFraction> elements = {silicon, carbon, titanium, copper};
    std::vector<ElementFraction> shuffled = {carbon, silicon, titanium, copper};

    /// create the material composition
    MaterialComposition elementsC(elements);
    MaterialComposition shuffledC(shuffled);
    // check if the sorting worked
    BOOST_CHECK_EQUAL(elementsC.size(), shuffledC.size());
    BOOST_CHECK_EQUAL(elementsC[0].first, shuffledC[0].first);
    BOOST_CHECK_EQUAL(elementsC[1].first, shuffledC[1].first);
    BOOST_CHECK_EQUAL(elementsC[2].first, shuffledC[2].first);
    BOOST_CHECK_EQUAL(elementsC[3].first, shuffledC[3].first);
    BOOST_CHECK_EQUAL(elementsC[0].second, shuffledC[0].second);
    BOOST_CHECK_EQUAL(elementsC[1].second, shuffledC[1].second);
    BOOST_CHECK_EQUAL(elementsC[2].second, shuffledC[2].second);
    BOOST_CHECK_EQUAL(elementsC[3].second, shuffledC[3].second);

    /// @todo the fraction test will only work when re-scaling is implemented
    /// check if the total fraction is one
    // float totalFraction = 0.;
    // for (auto& eFraction : elementsC)
    //  totalFraction += eFraction.fraction();
    // BOOST_CHECK_CLOSE(1., totalFraction, 2*elMaxTolerance);
  }

  BOOST_AUTO_TEST_CASE(Material_boolean_test)
  {

    Material vacuum;
    BOOST_CHECK_EQUAL(bool(vacuum), false);

    Material something(1., 2., 3., 4., 5);
    BOOST_CHECK_EQUAL(bool(something), true);
  }
}
}

// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Material Composition Tests
#include <boost/test/included/unit_test.hpp>
// clang-format on

#include <climits>

#include "Acts/Material/Material.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

namespace Test {

  // the maximum tolerance is half the accuracy
  float elMaxTolerance = 0.5 / float(UCHAR_MAX);

  namespace au = Acts::units;

  BOOST_AUTO_TEST_CASE(ElementFraction_construction_test)
  {
    // carbon parameters, atomic charge is Z
    unsigned char carbonZ = 12;
    // a fraction between 0 and 255
    unsigned char carbonFractionCH    = 46;
    float         carbonFractionFloat = float(46. / UCHAR_MAX);

    // test the carbon fraction
    std::array<unsigned char, 2> data = {{carbonZ, carbonFractionCH}};
    ElementFraction carbonEFC(data);
    ElementFraction carbonEFF((unsigned int)carbonZ, carbonFractionFloat);

    // check if you get the element and the fraction back
    BOOST_CHECK_EQUAL(12ul, carbonEFC.element());
    CHECK_CLOSE_REL(carbonFractionFloat, carbonEFC.fraction(), elMaxTolerance);
    BOOST_CHECK_EQUAL(12ul, carbonEFF.element());
    CHECK_CLOSE_REL(carbonFractionFloat, carbonEFF.fraction(), elMaxTolerance);
  }

  BOOST_AUTO_TEST_CASE(ElementFraction_movable_test)
  {
    // carbon parameters, atomic charge is Z
    unsigned char carbonZ = 12;
    // a fraction between 0 and 255
    float carbonFraction = float(46. / UCHAR_MAX);

    // test the carbon fraction
    ElementFraction carbon((unsigned int)carbonZ, carbonFraction);
    ElementFraction carbonMoved(std::move(carbon));

    BOOST_CHECK_EQUAL(12ul, carbonMoved.element());
    CHECK_CLOSE_REL(carbonFraction, carbonMoved.fraction(), elMaxTolerance);

    ElementFraction carbonMovedAssigned = std::move(carbonMoved);
    BOOST_CHECK_EQUAL(12ul, carbonMovedAssigned.element());
    CHECK_CLOSE_REL(
        carbonFraction, carbonMovedAssigned.fraction(), elMaxTolerance);
  }

  BOOST_AUTO_TEST_CASE(MaterialComposition_construction_test)
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
    BOOST_CHECK_EQUAL(elementsC.elements()[0].data()[0],
                      shuffledC.elements()[0].data()[0]);
    BOOST_CHECK_EQUAL(elementsC.elements()[1].data()[0],
                      shuffledC.elements()[1].data()[0]);
    BOOST_CHECK_EQUAL(elementsC.elements()[2].data()[0],
                      shuffledC.elements()[2].data()[0]);
    BOOST_CHECK_EQUAL(elementsC.elements()[3].data()[0],
                      shuffledC.elements()[3].data()[0]);
    BOOST_CHECK_EQUAL(elementsC.elements()[0].data()[1],
                      shuffledC.elements()[0].data()[1]);
    BOOST_CHECK_EQUAL(elementsC.elements()[1].data()[1],
                      shuffledC.elements()[1].data()[1]);
    BOOST_CHECK_EQUAL(elementsC.elements()[2].data()[1],
                      shuffledC.elements()[2].data()[1]);
    BOOST_CHECK_EQUAL(elementsC.elements()[3].data()[1],
                      shuffledC.elements()[3].data()[1]);
    /// or shortly
    BOOST_CHECK_EQUAL(elementsC, shuffledC);

    float totalFraction = 0.;
    for (auto& eFraction : elementsC.elements()) {
      totalFraction += eFraction.fraction();
    }
    // to better fit we need to implement some proper weight scaling
    BOOST_CHECK_LT(std::abs(1. - totalFraction),
                   elementsC.elements().size() * elMaxTolerance);
  }

  BOOST_AUTO_TEST_CASE(MaterialComposition_movable_test)
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

    /// create the material composition - elements are copied here
    MaterialComposition mComposition(elements);

    /// Move copy construction - object is moved here
    MaterialComposition mCompositionMoved(std::move(mComposition));
    BOOST_CHECK_EQUAL(mCompositionMoved.size(), elements.size());

    MaterialComposition mCompositionMovedAssigned
        = std::move(mCompositionMoved);
    BOOST_CHECK_EQUAL(mCompositionMovedAssigned.size(), elements.size());
  }
}
}
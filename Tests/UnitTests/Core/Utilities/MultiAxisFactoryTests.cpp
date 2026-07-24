// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/AxisFactory.hpp"
#include "Acts/Utilities/MultiAxisFactory.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <stdexcept>

using Acts::AxisFactory;
using Acts::AxisResolution;
using Acts::MultiAxisFactory;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(MultiAxisFactoryBasics) {
  using enum Acts::AxisBoundaryType;

  MultiAxisFactory maf({AxisFactory::Equidistant(Bound, 0., 1., 10),
                        AxisFactory::Variable(Open, {0., 1., 3.})});
  BOOST_CHECK_EQUAL(maf.size(), 2);
  BOOST_CHECK(!maf.isDeferred());
  BOOST_CHECK_EQUAL(maf.axisFactory(0).nBins(), 10);
  BOOST_CHECK_EQUAL(maf.axisFactory(1).nBins(), 2);
  BOOST_CHECK_THROW(maf.axisFactory(2), std::out_of_range);

  auto axes = maf.toAxes();
  BOOST_CHECK_EQUAL(axes.size(), 2);
  BOOST_CHECK(axes[0]->isEquidistant());
  BOOST_CHECK(axes[1]->isVariable());

  auto multiAxis = maf.toMultiAxis();
  BOOST_CHECK_EQUAL(multiAxis->getNAxes(), 2);
  BOOST_CHECK_EQUAL(multiAxis->getAxis(0).getNBins(), 10);
  BOOST_CHECK_EQUAL(multiAxis->getAxis(1).getNBins(), 2);

  // Empty construction is invalid
  BOOST_CHECK_THROW(MultiAxisFactory({}), std::invalid_argument);

  // Mixed with a deferred axis: full resolution is required for all axes
  MultiAxisFactory mixed({AxisFactory::Equidistant(Bound, 0., 1., 10),
                          AxisFactory::DeferredEquidistant(5)});
  BOOST_CHECK(mixed.isDeferred());
  BOOST_CHECK_THROW(mixed.toAxes(), std::domain_error);
}

BOOST_AUTO_TEST_CASE(MultiAxisFactoryDeferredResolution) {
  using enum Acts::AxisBoundaryType;
  using enum Acts::AxisDirection;

  MultiAxisFactory maf({AxisFactory::DeferredEquidistant(4, AxisRPhi),
                        AxisFactory::DeferredVariable({0., 0.5, 1.}, AxisZ)});
  BOOST_CHECK(maf.isDeferred());
  BOOST_CHECK_THROW(maf.toAxes(), std::domain_error);

  std::vector<AxisResolution> resolutions = {{-3., 3., Closed},
                                             {-10., 10., Bound}};
  auto axes = maf.toAxes(resolutions);
  BOOST_CHECK_EQUAL(axes.size(), 2);
  BOOST_CHECK_EQUAL(axes[0]->getBoundaryType(), Closed);
  BOOST_CHECK_EQUAL(axes[0]->getNBins(), 4);
  BOOST_CHECK(axes[0]->getDirection() == AxisRPhi);
  BOOST_CHECK_EQUAL(axes[1]->getBoundaryType(), Bound);
  CHECK_CLOSE_ABS(axes[1]->getBinEdges()[1], 0., 1e-15);
  BOOST_CHECK(axes[1]->getDirection() == AxisZ);

  // Directions are validated per axis
  std::vector<Acts::AxisDirection> good = {AxisRPhi, AxisZ};
  std::vector<Acts::AxisDirection> bad = {AxisZ, AxisRPhi};
  BOOST_CHECK_NO_THROW(maf.toAxes(resolutions, good));
  BOOST_CHECK_THROW(maf.toAxes(resolutions, bad), std::invalid_argument);

  // Sizes have to match
  std::vector<AxisResolution> tooFew = {{0., 1., Bound}};
  BOOST_CHECK_THROW(maf.toAxes(tooFew), std::invalid_argument);
  std::vector<Acts::AxisDirection> tooFewDirs = {AxisRPhi};
  BOOST_CHECK_THROW(maf.toAxes(resolutions, tooFewDirs), std::invalid_argument);

  auto multiAxis = maf.toMultiAxis(resolutions);
  BOOST_CHECK_EQUAL(multiAxis->getNAxes(), 2);
  BOOST_CHECK_EQUAL(multiAxis->getAxis(0).getBoundaryType(), Closed);
}

BOOST_AUTO_TEST_CASE(MultiAxisFactoryXDApi) {
  using enum Acts::AxisBoundaryType;

  Acts::MultiAxisFactory1D maf1D({AxisFactory::DeferredEquidistant(8)});
  BOOST_CHECK_EQUAL(maf1D.size(), 1);
  std::unique_ptr<Acts::IMultiAxis1D> ma1D =
      maf1D.toMultiAxis({AxisResolution{0., 4., Bound}});
  BOOST_CHECK_EQUAL(ma1D->getNAxes(), 1);
  BOOST_CHECK_EQUAL(ma1D->getAxis(0).getNBins(), 8);

  Acts::MultiAxisFactory2D maf2D({AxisFactory::Equidistant(Bound, 0., 1., 2),
                                  AxisFactory::Equidistant(Bound, 0., 1., 3)});
  std::unique_ptr<Acts::IMultiAxis2D> ma2D = maf2D.toMultiAxis();
  BOOST_CHECK_EQUAL(ma2D->getNAxes(), 2);
  BOOST_CHECK_EQUAL(ma2D->getNTotalBins(), 6);

  // The XD types remain usable through the runtime-dimension base
  const MultiAxisFactory& base = maf2D;
  BOOST_CHECK_EQUAL(base.toMultiAxis()->getNAxes(), 2);
}

BOOST_AUTO_TEST_CASE(MultiAxisFactoryEqualityAndStreams) {
  using enum Acts::AxisBoundaryType;

  MultiAxisFactory a({AxisFactory::DeferredEquidistant(4),
                      AxisFactory::DeferredEquidistant(5)});
  MultiAxisFactory b({AxisFactory::DeferredEquidistant(4),
                      AxisFactory::DeferredEquidistant(5)});
  MultiAxisFactory c({AxisFactory::DeferredEquidistant(4)});
  BOOST_CHECK(a == b);
  BOOST_CHECK(a != c);

  BOOST_CHECK_EQUAL(
      c.toString(),
      "MultiAxisFactory: 1 axes [AxisFactory: 4 bins, equidistant within "
      "deferred range]");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

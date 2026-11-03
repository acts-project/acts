// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// GCC 15 reports the axis metadata strings of the default `ResPlotTool::Config`
// as maybe-uninitialized while copying them into the binning map. The pragma
// has to cover the includes because the warning is reported at the inlined
// `std::string` destructor in `<string>`.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "ActsExamples/Validation/ResPlotTool.hpp"

#include <memory>
#include <numbers>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

using namespace Acts;
using namespace ActsExamples;

namespace {

/// The number of entries in the in-range bins of @p histogram, i.e. excluding
/// the under- and overflow bins.
double inRangeEntries(const ResPlotTool::Histogram1& histogram) {
  double sum = 0;
  for (int i = 0; i < histogram.histogram().axis(0).size(); ++i) {
    sum += histogram.binContent({i});
  }
  return sum;
}

/// The default config with the perigee parameter names filled in, which the
/// tool leaves to its caller.
ResPlotTool::Config makeConfig() {
  ResPlotTool::Config cfg;
  cfg.paramNames = {"d0", "z0", "phi", "theta", "qop", "t"};
  return cfg;
}

BoundTrackParameters makeParameters(
    const std::shared_ptr<const Surface>& surface, double phi) {
  BoundVector parameters = BoundVector::Zero();
  parameters[eBoundPhi] = phi;
  parameters[eBoundTheta] = std::numbers::pi / 2;
  parameters[eBoundQOverP] = 1. / (1. * UnitConstants::GeV);

  return BoundTrackParameters(surface, parameters,
                              BoundMatrix::Identity() * 1e-4,
                              ParticleHypothesis::pion());
}

}  // namespace

BOOST_AUTO_TEST_SUITE(ValidationResPlotTool)

// a track pointing along -x has truth and reco phi on opposite sides of the
// wrap, where the naive difference is almost a full turn
BOOST_AUTO_TEST_CASE(PhiResidualWrapsAround) {
  ResPlotTool tool(makeConfig(), Acts::Logging::INFO);

  const auto surface = Surface::makeShared<PerigeeSurface>(Vector3::Zero());
  const double delta = 0.00105;
  tool.fill(makeParameters(surface, std::numbers::pi - delta),
            makeParameters(surface, -std::numbers::pi + delta));

  // without the wrap the residual is `-2 pi + 2 delta` and lands in underflow
  const ResPlotTool::Histogram1& residual = tool.res().at("phi");
  BOOST_CHECK_EQUAL(inRangeEntries(residual), 1.);

  const int bin = residual.histogram().axis(0).index(2 * delta);
  BOOST_CHECK_EQUAL(residual.binContent({bin}), 1.);
}

// the same fill away from the wrap must be untouched by it
BOOST_AUTO_TEST_CASE(PhiResidualAwayFromTheWrap) {
  ResPlotTool tool(makeConfig(), Acts::Logging::INFO);

  const auto surface = Surface::makeShared<PerigeeSurface>(Vector3::Zero());
  const double delta = 0.00105;
  tool.fill(makeParameters(surface, 0.), makeParameters(surface, delta));

  const ResPlotTool::Histogram1& residual = tool.res().at("phi");
  const int bin = residual.histogram().axis(0).index(delta);
  BOOST_CHECK_EQUAL(residual.binContent({bin}), 1.);
}

BOOST_AUTO_TEST_SUITE_END()

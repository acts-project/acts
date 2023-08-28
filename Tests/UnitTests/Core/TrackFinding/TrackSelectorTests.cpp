// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/TrackFinding/TrackSelector.hpp"

using namespace Acts;

struct MockTrack {
  double m_theta;
  double m_phi;
  double m_pt;
  double m_loc0;
  double m_loc1;
  double m_time;
  double m_nMeasurements;

  double theta() const { return m_theta; }
  double phi() const { return m_phi; }
  double transverseMomentum() const { return m_pt; }
  double loc0() const { return m_loc0; }
  double loc1() const { return m_loc1; }
  double time() const { return m_time; }
  double nMeasurements() const { return m_nMeasurements; }
};

BOOST_AUTO_TEST_SUITE(TrackSelectorTests)

BOOST_AUTO_TEST_CASE(TestSingleBinCase) {
  TrackSelector::Config cfgBase;

  auto thetaFromEta = [](double eta) { return 2 * std::atan(std::exp(-eta)); };

  MockTrack baseTrack{};
  baseTrack.m_theta = thetaFromEta(0.5);
  baseTrack.m_phi = 0.5;
  baseTrack.m_pt = 0.5;
  baseTrack.m_loc0 = 0.5;
  baseTrack.m_loc1 = 0.5;
  baseTrack.m_time = 0.5;
  baseTrack.m_nMeasurements = 0.5;

  {
    TrackSelector selector{cfgBase};
    // should select anything
    BOOST_CHECK(selector.isValidTrack(baseTrack));
  }

  auto check = [&](const auto& var, const auto& minPtr, const auto& maxPtr,
                   const auto& propPtr) {
    BOOST_TEST_INFO_SCOPE("Testing " << var);
    MockTrack track = baseTrack;

    auto cfgMinOnly = cfgBase;
    auto cfgMaxOnly = cfgBase;
    auto cfgMinMax = cfgBase;

    cfgMinOnly.*minPtr = {-1};
    cfgMinMax.*minPtr = {-1};
    cfgMaxOnly.*maxPtr = {1};
    cfgMinMax.*maxPtr = {1};

    TrackSelector minOnly{cfgMinOnly};
    TrackSelector maxOnly{cfgMaxOnly};
    TrackSelector minMax{cfgMinMax};

    BOOST_CHECK(minOnly.isValidTrack(track));
    BOOST_CHECK(maxOnly.isValidTrack(track));
    BOOST_CHECK(minMax.isValidTrack(track));

    // push track outside of minimum
    track.*propPtr = -1.1;

    BOOST_CHECK(!minOnly.isValidTrack(track));
    BOOST_CHECK(maxOnly.isValidTrack(track));
    BOOST_CHECK(!minMax.isValidTrack(track));

    // push track outside of maximum
    track.*propPtr = 1.1;

    BOOST_CHECK(minOnly.isValidTrack(track));
    BOOST_CHECK(!maxOnly.isValidTrack(track));
    BOOST_CHECK(!minMax.isValidTrack(track));
  };

  check("loc0", &TrackSelector::Config::loc0Min,
        &TrackSelector::Config::loc0Max, &MockTrack::m_loc0);

  check("loc1", &TrackSelector::Config::loc1Min,
        &TrackSelector::Config::loc1Max, &MockTrack::m_loc1);

  check("phi", &TrackSelector::Config::phiMin, &TrackSelector::Config::phiMax,
        &MockTrack::m_phi);

  check("time", &TrackSelector::Config::timeMin,
        &TrackSelector::Config::timeMax, &MockTrack::m_time);

  {
    BOOST_TEST_INFO_SCOPE("pt min");
    auto cfg = cfgBase;
    cfg.ptMin = {0.2};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_pt = 0.1;
    BOOST_CHECK(!selector.isValidTrack(track));
  }

  {
    BOOST_TEST_INFO_SCOPE("pt max");
    auto cfg = cfgBase;
    cfg.ptMax = {1.0};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_pt = 1.1;
    BOOST_CHECK(!selector.isValidTrack(track));
  }

  {
    BOOST_TEST_INFO_SCOPE("pt min max");
    auto cfg = cfgBase;
    cfg.ptMin = {0.2};
    cfg.ptMax = {1.0};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_pt = 0.1;
    BOOST_CHECK(!selector.isValidTrack(track));
    track.m_pt = 1.1;
    BOOST_CHECK(!selector.isValidTrack(track));
  }

  {
    BOOST_TEST_INFO_SCOPE("eta min");
    auto cfg = cfgBase;
    cfg.etaMin = {-1.0};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_theta = thetaFromEta(-1.1);
    BOOST_CHECK(!selector.isValidTrack(track));
  }

  {
    BOOST_TEST_INFO_SCOPE("eta max");
    auto cfg = cfgBase;
    cfg.etaMax = {1.0};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_theta = thetaFromEta(1.1);
    BOOST_CHECK(!selector.isValidTrack(track));
  }

  {
    BOOST_TEST_INFO_SCOPE("eta min max");
    auto cfg = cfgBase;
    cfg.etaMin = {-1.0};
    cfg.etaMax = {1.0};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_theta = thetaFromEta(-1.1);
    BOOST_CHECK(!selector.isValidTrack(track));
    track.m_theta = thetaFromEta(1.1);
    BOOST_CHECK(!selector.isValidTrack(track));
  }

  {
    BOOST_TEST_INFO_SCOPE("abs eta min");
    auto cfg = cfgBase;
    cfg.absEtaMin = {0.2};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_theta = thetaFromEta(-0.5);
    BOOST_CHECK(selector.isValidTrack(track));

    track.m_theta = thetaFromEta(0.1);
    BOOST_CHECK(!selector.isValidTrack(track));
    track.m_theta = thetaFromEta(-0.1);
    BOOST_CHECK(!selector.isValidTrack(track));
  }

  {
    BOOST_TEST_INFO_SCOPE("abs eta max");
    auto cfg = cfgBase;
    cfg.absEtaMax = {1.0};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_theta = thetaFromEta(-0.5);
    BOOST_CHECK(selector.isValidTrack(track));

    track.m_theta = thetaFromEta(1.1);
    BOOST_CHECK(!selector.isValidTrack(track));
    track.m_theta = thetaFromEta(-1.1);
    BOOST_CHECK(!selector.isValidTrack(track));
  }

  {
    BOOST_TEST_INFO_SCOPE("abs eta min max");
    auto cfg = cfgBase;
    cfg.absEtaMin = {0.2};
    cfg.absEtaMax = {1.0};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_theta = thetaFromEta(-0.5);
    BOOST_CHECK(selector.isValidTrack(track));

    track.m_theta = thetaFromEta(0.1);
    BOOST_CHECK(!selector.isValidTrack(track));
    track.m_theta = thetaFromEta(-0.1);
    BOOST_CHECK(!selector.isValidTrack(track));

    track.m_theta = thetaFromEta(1.1);
    BOOST_CHECK(!selector.isValidTrack(track));
    track.m_theta = thetaFromEta(-1.1);
    BOOST_CHECK(!selector.isValidTrack(track));
  }

  {
    BOOST_TEST_INFO_SCOPE("nMeas min");
    auto cfg = cfgBase;
    cfg.minMeasurements = {1};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    track.m_nMeasurements = {2};
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_nMeasurements = {1};
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_nMeasurements = {0};
    BOOST_CHECK(!selector.isValidTrack(track));
  }
}

BOOST_AUTO_TEST_SUITE_END()

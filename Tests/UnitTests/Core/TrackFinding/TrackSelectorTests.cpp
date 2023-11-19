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

#include <limits>

using namespace Acts;
namespace bdata = boost::unit_test::data;

struct MockTrack {
  double m_theta;
  double m_phi;
  double m_pt;
  double m_loc0;
  double m_loc1;
  double m_time;
  double m_nMeasurements;

  bool hasReferenceSurface() const { return true; }
  double theta() const { return m_theta; }
  double phi() const { return m_phi; }
  double transverseMomentum() const { return m_pt; }
  double loc0() const { return m_loc0; }
  double loc1() const { return m_loc1; }
  double time() const { return m_time; }
  double nMeasurements() const { return m_nMeasurements; }
};

double thetaFromEta(double eta) {
  return 2 * std::atan(std::exp(-eta));
}

BOOST_AUTO_TEST_SUITE(TrackSelectorTests)

std::vector<double> etaValues{-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5,
                              -1.0, -0.5, 0.0,  0.5,  1.0,  1.5,  2.0,  2.5,
                              3.0,  3.5,  4.0,  4.5,  5.0,  1.0};

BOOST_DATA_TEST_CASE(TestSingleBinCase, bdata::make(etaValues), eta) {
  TrackSelector::EtaBinnedConfig cfgBase;

  MockTrack baseTrack{};
  baseTrack.m_theta = thetaFromEta(eta);
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

    cfgMinOnly.cutSets.at(0).*minPtr = -1;
    cfgMinMax.cutSets.at(0).*minPtr = -1;
    cfgMaxOnly.cutSets.at(0).*maxPtr = 1;
    cfgMinMax.cutSets.at(0).*maxPtr = 1;

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
    cfg.cutSets.at(0).ptMin = {0.2};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_pt = 0.1;
    BOOST_CHECK(!selector.isValidTrack(track));
  }

  {
    BOOST_TEST_INFO_SCOPE("pt max");
    auto cfg = cfgBase;
    cfg.cutSets.at(0).ptMax = {1.0};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_pt = 1.1;
    BOOST_CHECK(!selector.isValidTrack(track));
  }

  {
    BOOST_TEST_INFO_SCOPE("pt min max");
    auto cfg = cfgBase;
    cfg.cutSets.at(0).ptMin = {0.2};
    cfg.cutSets.at(0).ptMax = {1.0};
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
    cfg.cutSets.at(0).etaMin = {-1.0};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    track.m_theta = thetaFromEta(0.5);
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_theta = thetaFromEta(-1.1);
    BOOST_CHECK(!selector.isValidTrack(track));
  }

  {
    BOOST_TEST_INFO_SCOPE("eta max");
    auto cfg = cfgBase;
    cfg.cutSets.at(0).etaMax = {1.0};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    track.m_theta = thetaFromEta(0.5);
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_theta = thetaFromEta(1.1);
    BOOST_CHECK(!selector.isValidTrack(track));
  }

  {
    BOOST_TEST_INFO_SCOPE("eta min max");
    auto cfg = cfgBase;
    cfg.cutSets.at(0).etaMin = {-1.0};
    cfg.cutSets.at(0).etaMax = {1.0};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    track.m_theta = thetaFromEta(0.5);
    BOOST_CHECK(selector.isValidTrack(track));
    track.m_theta = thetaFromEta(-1.1);
    BOOST_CHECK(!selector.isValidTrack(track));
    track.m_theta = thetaFromEta(1.1);
    BOOST_CHECK(!selector.isValidTrack(track));
  }

  {
    BOOST_TEST_INFO_SCOPE("abs eta min");
    auto cfg = cfgBase;
    cfg.cutSets.at(0).absEtaMin = {0.2};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    track.m_theta = thetaFromEta(0.5);
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
    cfg.cutSets.at(0).absEtaMax = {1.0};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    track.m_theta = thetaFromEta(0.5);
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
    cfg.cutSets.at(0).absEtaMin = {0.2};
    cfg.cutSets.at(0).absEtaMax = {1.0};
    TrackSelector selector{cfg};
    MockTrack track = baseTrack;
    track.m_theta = thetaFromEta(0.5);
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
    cfg.cutSets.at(0).minMeasurements = {1};
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

BOOST_AUTO_TEST_CASE(TestSingleBinEtaCutByBinEdge) {
  TrackSelector selector{TrackSelector::EtaBinnedConfig(1.0).addCuts(2.0)};

  BOOST_TEST_INFO_SCOPE(selector.config());

  MockTrack track{};
  track.m_theta = thetaFromEta(0.0);
  BOOST_CHECK(!selector.isValidTrack(track));

  track.m_theta = thetaFromEta(0.5);
  BOOST_CHECK(!selector.isValidTrack(track));

  // cannot easily check on-edge behavior because of floating point arithmetic
  // (it won't be exactly 1.0 in selector)
  track.m_theta = thetaFromEta(1.01);
  BOOST_CHECK(selector.isValidTrack(track));

  track.m_theta = thetaFromEta(1.5);
  BOOST_CHECK(selector.isValidTrack(track));

  track.m_theta = thetaFromEta(2.0);
  BOOST_CHECK(!selector.isValidTrack(track));
}

BOOST_AUTO_TEST_CASE(TestMultiBinCuts) {
  MockTrack baseTrack{};
  baseTrack.m_theta = thetaFromEta(1.0);
  baseTrack.m_phi = 0.5;
  baseTrack.m_pt = 0.5;
  baseTrack.m_loc0 = 0.5;
  baseTrack.m_loc1 = 0.5;
  baseTrack.m_time = 0.5;
  baseTrack.m_nMeasurements = 0.5;

  using Config = TrackSelector::Config;

  using factory_ptr_t = Config& (Config::*)(double, double);
  using prop_ptr_t = double MockTrack::*;

  auto check = [&](const char* name, const factory_ptr_t& factory,
                   const prop_ptr_t& prop) {
    BOOST_TEST_CONTEXT(name) {
      auto cfg = TrackSelector::EtaBinnedConfig{0.0};

      cfg.addCuts(2.0, [&](auto& c) { (c.*factory)(-1.0, 1.0); })
          .addCuts([&](auto& c) { (c.*factory)(-2.0, 2.0); });

      TrackSelector selector{cfg};

      BOOST_TEST_INFO_SCOPE(cfg);

      {
        // exactly at zero
        MockTrack track = baseTrack;
        track.m_theta = thetaFromEta(0.0);

        BOOST_CHECK(selector.isValidTrack(track));

        track.*prop = -1.1;
        BOOST_CHECK(!selector.isValidTrack(track));

        track.*prop = 1.1;
        BOOST_CHECK(!selector.isValidTrack(track));
      }

      {
        // first bin
        MockTrack track = baseTrack;
        track.m_theta = thetaFromEta(1.0);

        BOOST_CHECK(selector.isValidTrack(track));

        track.*prop = -1.1;
        BOOST_CHECK(!selector.isValidTrack(track));

        track.*prop = 1.1;
        BOOST_CHECK(!selector.isValidTrack(track));
      }

      {
        // first bin edge
        MockTrack track = baseTrack;
        track.m_theta =
            thetaFromEta(2.0 - std::numeric_limits<double>::epsilon());

        BOOST_CHECK(selector.isValidTrack(track));

        track.*prop = -1.1;
        BOOST_CHECK(!selector.isValidTrack(track));

        track.*prop = 1.1;
        BOOST_CHECK(!selector.isValidTrack(track));
      }

      {
        // second bin lower edge
        MockTrack track = baseTrack;
        track.m_theta = thetaFromEta(2.0);

        BOOST_CHECK(selector.isValidTrack(track));

        track.*prop = -1.1;
        BOOST_CHECK(selector.isValidTrack(track));

        track.*prop = 1.1;
        BOOST_CHECK(selector.isValidTrack(track));

        track.*prop = -2.1;
        BOOST_CHECK(!selector.isValidTrack(track));

        track.*prop = 2.1;
        BOOST_CHECK(!selector.isValidTrack(track));
      }

      {
        // second bin
        MockTrack track = baseTrack;
        track.m_theta = thetaFromEta(666.0);

        track.*prop = -1.1;
        BOOST_CHECK(selector.isValidTrack(track));

        track.*prop = 1.1;
        BOOST_CHECK(selector.isValidTrack(track));

        track.*prop = -2.1;
        BOOST_CHECK(!selector.isValidTrack(track));

        track.*prop = 2.1;
        BOOST_CHECK(!selector.isValidTrack(track));
      }
    }
  };

  check("loc0", &Config::loc0, &MockTrack::m_loc0);
  check("loc1", &Config::loc1, &MockTrack::m_loc1);
  check("time", &Config::time, &MockTrack::m_time);
  check("phi", &Config::phi, &MockTrack::m_phi);
  check("pt", &Config::pt, &MockTrack::m_pt);
}

BOOST_AUTO_TEST_CASE(TestBinSelection) {
  using EtaBinnedConfig = TrackSelector::EtaBinnedConfig;
  constexpr double inf = std::numeric_limits<double>::infinity();

  {
    EtaBinnedConfig cfg{std::vector<double>{0, inf}};
    for (int i = -1; i <= 1; i = i + 2) {
      BOOST_CHECK_EQUAL(cfg.binIndex(i * 0.0), 0);
      BOOST_CHECK_EQUAL(cfg.binIndex(i * 1.0), 0);
      BOOST_CHECK_EQUAL(cfg.binIndex(i * 2.0), 0);
      BOOST_CHECK_EQUAL(cfg.binIndex(i * 3.0), 0);
      BOOST_CHECK_EQUAL(cfg.binIndex(i * 10.0), 0);
      BOOST_CHECK_EQUAL(cfg.binIndex(i * 1000.0), 0);
    }
  }

  {
    EtaBinnedConfig cfg{std::vector<double>{0, 0.5, 1.5, 2.5, 3.0, inf}};
    for (int i = -1; i <= 1; i = i + 2) {
      BOOST_CHECK_EQUAL(cfg.binIndex(i * 0.0), 0);
      BOOST_CHECK_EQUAL(cfg.binIndex(i * 1.0), 1);
      BOOST_CHECK_EQUAL(cfg.binIndex(i * 2.0), 2);
      BOOST_CHECK_EQUAL(cfg.binIndex(i * 3.0), 4);
      BOOST_CHECK_EQUAL(cfg.binIndex(i * 10.0), 4);
      BOOST_CHECK_EQUAL(cfg.binIndex(i * 1000.0), 4);
    }
  }

  {
    EtaBinnedConfig cfg{std::vector<double>{0, 1, 2}};
    for (int i = -1; i <= 1; i = i + 2) {
      BOOST_CHECK_EQUAL(cfg.binIndex(i * 0.0), 0);
      BOOST_CHECK_EQUAL(cfg.binIndex(i * 1.0), 1);
      BOOST_CHECK_EQUAL(
          cfg.binIndex(i * (2.0 - std::numeric_limits<double>::epsilon())), 1);
      BOOST_CHECK_THROW(cfg.binIndex(i * 2.0), std::invalid_argument);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestConstructor) {
  // valid multi bin construction
  {
    TrackSelector::EtaBinnedConfig cfg{std::vector<double>{0, 1, 4}};
    cfg.cutSets.at(0).ptMin = 0.9;
    cfg.cutSets.at(1).ptMin = 0.4;
    TrackSelector{cfg};
  }

  {
    // Track selector config with 2 eta bins
    TrackSelector::EtaBinnedConfig cfg{std::vector<double>{0, 1, 4}};

    // just the right amount!
    TrackSelector{cfg};

    // not enough cut cets
    cfg.cutSets.resize(1);
    BOOST_CHECK_THROW(TrackSelector{cfg}, std::invalid_argument);

    // too many cut sets
    cfg.cutSets.resize(3);
    BOOST_CHECK_THROW(TrackSelector{cfg}, std::invalid_argument);
  }

  // Constructor straight from cut config
  TrackSelector::Config cuts;
  TrackSelector selector{cuts};
  BOOST_CHECK_EQUAL(selector.config().absEtaEdges.size(), 2);

  {
    // Invalid sequence of chained construction
    auto cfg = TrackSelector::EtaBinnedConfig(0);
    cfg.addCuts(2.0, [](auto& c) { c.loc0(-1.0, 1.0); });
    BOOST_CHECK_THROW(cfg.addCuts(1.0), std::invalid_argument);
    BOOST_CHECK_THROW(TrackSelector::EtaBinnedConfig(0).addCuts(-2.0),
                      std::invalid_argument);
  }

  {
    auto cfg = TrackSelector::EtaBinnedConfig(1.0);

    cfg.addCuts(2.0, [](auto& c) { c.loc0(-1.0, 1.0); });
    BOOST_CHECK_EQUAL(cfg.nEtaBins(), 1);
    BOOST_CHECK_EQUAL(cfg.getCuts(1.5).loc0Min, -1.0);
    BOOST_CHECK_EQUAL(cfg.getCuts(1.5).loc0Max, 1.0);

    cfg.addCuts(3.0, [](auto& c) { c.loc0(-2.0, 2.0); });
    BOOST_CHECK_EQUAL(cfg.nEtaBins(), 2);
    BOOST_CHECK_EQUAL(cfg.getCuts(2.5).loc0Min, -2.0);
    BOOST_CHECK_EQUAL(cfg.getCuts(2.5).loc0Max, 2.0);

    cfg.addCuts(4.0, [](auto& c) { c.loc0(-3.0, 3.0); });
    BOOST_CHECK_EQUAL(cfg.nEtaBins(), 3);
    BOOST_CHECK_EQUAL(cfg.getCuts(3.5).loc0Min, -3.0);
    BOOST_CHECK_EQUAL(cfg.getCuts(3.5).loc0Max, 3.0);
  }
}

BOOST_AUTO_TEST_SUITE_END()

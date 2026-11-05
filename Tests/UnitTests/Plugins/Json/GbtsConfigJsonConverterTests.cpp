// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "ActsPlugins/Json/GbtsConfigJsonConverter.hpp"

#include <nlohmann/json.hpp>

using Acts::Experimental::GbtsLayerConfig;

BOOST_AUTO_TEST_SUITE(GbtsConfigJsonConverter)

BOOST_AUTO_TEST_CASE(LayerDepthIsOptional) {
  nlohmann::json j = {{"id", 8},
                      {"type", "barrel"},
                      {"technology", "pixel"},
                      {"surfaces", nlohmann::json::array()}};

  // a layer that does not name a depth is out of the innermost-layer cuts
  BOOST_CHECK_EQUAL(j.get<GbtsLayerConfig>().depth, -1);

  j["depth"] = 2;
  BOOST_CHECK_EQUAL(j.get<GbtsLayerConfig>().depth, 2);
}

BOOST_AUTO_TEST_CASE(LayerDepthRoundTrip) {
  GbtsLayerConfig layer;
  layer.id = 8;

  // left out when unset, so that a map without depths stays as it was
  nlohmann::json j = layer;
  BOOST_CHECK(!j.contains("depth"));
  BOOST_CHECK_EQUAL(j.get<GbtsLayerConfig>().depth, -1);

  layer.depth = 3;
  j = layer;
  BOOST_CHECK_EQUAL(j.get<GbtsLayerConfig>().depth, 3);
}

BOOST_AUTO_TEST_SUITE_END()

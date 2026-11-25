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

BOOST_AUTO_TEST_CASE(LayerDepthIsRequired) {
  nlohmann::json j = {{"id", 8},
                      {"type", "barrel"},
                      {"technology", "pixel"},
                      {"surfaces", nlohmann::json::array()}};

  // every layer must name its depth, -1 for outside the innermost-layer cuts
  BOOST_CHECK_THROW(j.get<GbtsLayerConfig>(), nlohmann::json::out_of_range);

  j["depth"] = -1;
  BOOST_CHECK_EQUAL(j.get<GbtsLayerConfig>().depth, -1);

  j["depth"] = 2;
  BOOST_CHECK_EQUAL(j.get<GbtsLayerConfig>().depth, 2);
}

BOOST_AUTO_TEST_CASE(LayerDepthRoundTrip) {
  GbtsLayerConfig layer;
  layer.id = 8;

  // always written, so a written map names the depth of every layer
  nlohmann::json j = layer;
  BOOST_CHECK_EQUAL(j.at("depth").get<std::int32_t>(), -1);
  BOOST_CHECK_EQUAL(j.get<GbtsLayerConfig>().depth, -1);

  layer.depth = 3;
  j = layer;
  BOOST_CHECK_EQUAL(j.get<GbtsLayerConfig>().depth, 3);
}

BOOST_AUTO_TEST_SUITE_END()

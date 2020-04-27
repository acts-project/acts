// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <ios>
#include <iostream>
#include <stdexcept>

#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Plugins/Json/JsonGeometryConverter.hpp"

using json = nlohmann::json;

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_CASE(Json_conversion) {
  std::stringstream ifj;

  Acts::JsonGeometryConverter::Config cfg;

  ifj << "{";
  ifj << "    \"volumes\": {";
  ifj << "        \"2\": {";
  ifj << "            \"Geoid\": \"[   2 |   0 |   2 |   0 |    0 ]\",";
  ifj << "            \"Name\": \"\",";
  ifj << "            \"layers\": {";
  ifj << "                \"2\": {";
  ifj << "                    \"Geoid\": \"[   2 |   0 |   2 |   0 |    0 ]\",";
  ifj << "                    \"representing\": {";
  ifj << "                        \"bin0\": [";
  ifj << "                            \"binR\",";
  ifj << "                            \"open\",";
  ifj << "                            2,";
  ifj << "                            [";
  ifj << "                                28.0,";
  ifj << "                                186.0";
  ifj << "                            ]";
  ifj << "                        ],";
  ifj << "                        \"data\": [";
  ifj << "                            [";
  ifj << "                                [";
  ifj << "                                    116.40576171875,";
  ifj << "                                    220.39578247070313,";
  ifj << "                                    7.115414142608643,";
  ifj << "                                    1.42858931477799e+22,";
  ifj << "                                    0.003504054620862007,";
  ifj << "                                    1.0";
  ifj << "                                ],";
  ifj << "                                [";
  ifj << "                                    431.8612976074219,";
  ifj << "                                    1103.45947265625,";
  ifj << "                                    11.447914123535156,";
  ifj << "                                    4.5373971722490746e+21,";
  ifj << "                                    0.0007596652721986175,";
  ifj << "                                    1.0";
  ifj << "                                ]";
  ifj << "                            ]";
  ifj << "                        ],";
  ifj << "                        \"mapMaterial\": true,";
  ifj << "                        \"type\": \"binned\"";
  ifj << "                    }";
  ifj << "                },";
  ifj << "                \"4\": {";
  ifj << "                    \"Geoid\": \"[   2 |   0 |   4 |   0 |    0 ]\",";
  ifj << "                    \"representing\": {";
  ifj << "                        \"bin0\": [";
  ifj << "                            \"binR\",";
  ifj << "                            \"open\",";
  ifj << "                            2,";
  ifj << "                            [";
  ifj << "                                28.0,";
  ifj << "                                186.0";
  ifj << "                            ]";
  ifj << "                        ],";
  ifj << "                        \"data\": [";
  ifj << "                            [";
  ifj << "                                [";
  ifj << "                                    91.93806457519531,";
  ifj << "                                    170.03167724609375,";
  ifj << "                                    4.947852611541748,";
  ifj << "                                    1.8783701146029598e+22,";
  ifj << "                                    0.004517487715929747,";
  ifj << "                                    1.0";
  ifj << "                                ],";
  ifj << "                                [";
  ifj << "                                    566.8548583984375,";
  ifj << "                                    1606.4842529296875,";
  ifj << "                                    10.814791679382324,";
  ifj << "                                    5.322017395554306e+21,";
  ifj << "                                    0.0005397787317633629,";
  ifj << "                                    1.0";
  ifj << "                                ]";
  ifj << "                            ]";
  ifj << "                        ],";
  ifj << "                        \"mapMaterial\": true,";
  ifj << "                        \"type\": \"binned\"";
  ifj << "                    }";
  ifj << "                }";
  ifj << "            }";
  ifj << "        }";
  ifj << "    }";
  ifj << "}";

  json jin;
  ifj >> jin;

  Acts::JsonGeometryConverter jmConverter(cfg);

  auto Map = jmConverter.jsonToMaterialMaps(jin);
  auto jint = jmConverter.materialMapsToJson(Map);

  BOOST_CHECK_EQUAL(jin, jint);
}

}  // namespace Test
}  // namespace Acts

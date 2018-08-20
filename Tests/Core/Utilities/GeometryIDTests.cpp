// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Boost include(s)
#define BOOST_TEST_MODULE GeometryID Tests
#include <boost/test/included/unit_test.hpp>
#include "Acts/Utilities/GeometryID.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {
namespace Test {

  geo_id_value volume_mask    = GeometryID::volume_mask;
  geo_id_value boundary_mask  = GeometryID::boundary_mask;
  geo_id_value layer_mask     = GeometryID::layer_mask;
  geo_id_value approach_mask  = GeometryID::approach_mask;
  geo_id_value sensitive_mask = GeometryID::sensitive_mask;
  geo_id_value channel_mask   = GeometryID::channel_mask;

  /// test of the geometry ID creation and consistency of the ranges
  BOOST_AUTO_TEST_CASE(GeometryID_test)
  {
    // create the volume shift and range
    geo_id_value volume_range = 64 - ACTS_BIT_SHIFT(volume_mask);

    geo_id_value boundary_range
        = ACTS_BIT_SHIFT(volume_mask) - ACTS_BIT_SHIFT(boundary_mask);

    geo_id_value layer_range
        = ACTS_BIT_SHIFT(boundary_mask) - ACTS_BIT_SHIFT(layer_mask);

    geo_id_value approach_range
        = ACTS_BIT_SHIFT(layer_mask) - ACTS_BIT_SHIFT(approach_mask);

    geo_id_value sensitive_range
        = ACTS_BIT_SHIFT(approach_mask) - ACTS_BIT_SHIFT(sensitive_mask);

    geo_id_value channel_range = ACTS_BIT_SHIFT(sensitive_mask)
        - ACTS_BIT_SHIFT(GeometryID::channel_mask);

    /// prepare all the masks and shifts
    std::vector<std::pair<geo_id_value, geo_id_value>> masks_range
        = {{volume_mask, volume_range},
           {boundary_mask, boundary_range},
           {layer_mask, layer_range},
           {approach_mask, approach_range},
           {sensitive_mask, sensitive_range},
           {channel_mask, channel_range}};

    for (auto msr : masks_range) {

      auto mask  = msr.first;
      auto range = msr.second;
      /// test range by [0, 1, 2^range-1]
      std::vector<geo_id_value> range_values
          = {0, 1, (geo_id_value(1) << range) - 1};
      for (auto& idv : range_values) {
        // create the geometry ID
        GeometryID geoID(idv, mask);
        // encode - decode test
        BOOST_CHECK_EQUAL(idv,
                          (ACTS_BIT_DECODE(ACTS_BIT_ENCODE(idv, mask), mask)));
        // geo id decoding
        BOOST_CHECK_EQUAL(idv, geoID.value(mask));
      }
    }
  }

  /// test the full encoding / decoding chain
  BOOST_AUTO_TEST_CASE(FullGeometryID_test)
  {
    // decode the IDs
    GeometryID volumeID(1, GeometryID::volume_mask);
    GeometryID boundaryID(2, GeometryID::boundary_mask);
    GeometryID layerID(3, GeometryID::layer_mask);
    GeometryID approachID(4, GeometryID::approach_mask);
    GeometryID sensitiveID(5, GeometryID::sensitive_mask);
    GeometryID channelID(6, GeometryID::channel_mask);

    // now check the validity before adding
    BOOST_CHECK_EQUAL(1lu, volumeID.value(GeometryID::volume_mask));
    BOOST_CHECK_EQUAL(2lu, boundaryID.value(GeometryID::boundary_mask));
    BOOST_CHECK_EQUAL(3lu, layerID.value(GeometryID::layer_mask));
    BOOST_CHECK_EQUAL(4lu, approachID.value(GeometryID::approach_mask));
    BOOST_CHECK_EQUAL(5lu, sensitiveID.value(GeometryID::sensitive_mask));
    BOOST_CHECK_EQUAL(6lu, channelID.value(GeometryID::channel_mask));

    // now create a compound ones
    GeometryID compoundID_dconst;
    compoundID_dconst += volumeID;
    GeometryID compoundID_cconst(volumeID);
    GeometryID compoundID_assign = volumeID;

    std::vector<GeometryID> compoundIDs
        = {compoundID_dconst, compoundID_cconst, compoundID_assign};

    /// check the validity after assigning/copying/constructing
    BOOST_CHECK_EQUAL(1lu, compoundID_dconst.value(GeometryID::volume_mask));
    BOOST_CHECK_EQUAL(1lu, compoundID_cconst.value(GeometryID::volume_mask));
    BOOST_CHECK_EQUAL(1lu, compoundID_assign.value(GeometryID::volume_mask));

    for (auto& cid : compoundIDs) {
      // add the sub IDs
      cid += boundaryID;
      cid += layerID;
      cid += approachID;
      cid += sensitiveID;
      cid += channelID;
      // now check the cid
      BOOST_CHECK_EQUAL(1lu, cid.value(GeometryID::volume_mask));
      BOOST_CHECK_EQUAL(2lu, cid.value(GeometryID::boundary_mask));
      BOOST_CHECK_EQUAL(3lu, cid.value(GeometryID::layer_mask));
      BOOST_CHECK_EQUAL(4lu, cid.value(GeometryID::approach_mask));
      BOOST_CHECK_EQUAL(5lu, cid.value(GeometryID::sensitive_mask));
      BOOST_CHECK_EQUAL(6lu, cid.value(GeometryID::channel_mask));
    }
  }

}  // namespace Test
}  // namespace Acts

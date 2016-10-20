// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Boost include(s)
#define BOOST_TEST_MODULE GeometryID Tests
#include "ACTS/Utilities/GeometryID.hpp"
#include <boost/test/included/unit_test.hpp>

namespace Acts {
namespace Test {

  /// prepare all the masks and shifts
  std::vector<std::vector<geo_id_value>> masks_shifts_range
      = {{GeometryID::volume_mask,
          GeometryID::volume_shift,
          GeometryID::volume_range},
         {GeometryID::boundary_mask,
          GeometryID::boundary_shift,
          GeometryID::boundary_range},
         {GeometryID::layer_mask,
          GeometryID::layer_shift,
          GeometryID::layer_range},
         {GeometryID::approach_mask,
          GeometryID::approach_shift,
          GeometryID::approach_range},
         {GeometryID::sensitive_mask,
          GeometryID::sensitive_shift,
          GeometryID::sensitive_range},
         {GeometryID::channel_mask,
          GeometryID::channel_shift,
          GeometryID::channel_range}};

  /// test of the geometry ID creation and consistency of the ranges
  BOOST_AUTO_TEST_CASE(GeometryID_test)
  {
    for (auto msr : masks_shifts_range) {
      auto mask  = msr[0];
      auto shift = msr[1];
      auto range = msr[2];
      /// test the full range of ids
      for (geo_id_value idv = 0; idv < pow(2, range); ++idv) {
        /// create the geometry ID
        GeometryID geoID(idv, shift);
        BOOST_CHECK_EQUAL(idv, geoID.value(mask, shift));
      }
    }
  }

  /// test the full encoding / decoding chain
  BOOST_AUTO_TEST_CASE(FullGeometryID_test)
  {
    // decode the IDs
    GeometryID volumeID(1, GeometryID::volume_shift);
    GeometryID boundaryID(2, GeometryID::boundary_shift);
    GeometryID layerID(3, GeometryID::layer_shift);
    GeometryID approachID(4, GeometryID::approach_shift);
    GeometryID sensitiveID(5, GeometryID::sensitive_shift);
    GeometryID channelID(6, GeometryID::channel_shift);
    // now create a compound ones
    GeometryID compoundID_dconst;
    compoundID_dconst += volumeID;
    GeometryID compoundID_cconst(volumeID);
    GeometryID compoundID_assign = volumeID;
    //
    std::vector<GeometryID> compoundIDs
        = {compoundID_dconst, compoundID_cconst, compoundID_assign};
    for (auto& cid : compoundIDs) {
      // add the sub IDs
      cid += boundaryID;
      cid += layerID;
      cid += approachID;
      cid += sensitiveID;
      cid += channelID;
      // now check the cid
      BOOST_CHECK_EQUAL(
          1lu, cid.value(GeometryID::volume_mask, GeometryID::volume_shift));
      BOOST_CHECK_EQUAL(
          2lu,
          cid.value(GeometryID::boundary_mask, GeometryID::boundary_shift));
      BOOST_CHECK_EQUAL(
          3lu, cid.value(GeometryID::layer_mask, GeometryID::layer_shift));
      BOOST_CHECK_EQUAL(
          4lu,
          cid.value(GeometryID::approach_mask, GeometryID::approach_shift));
      BOOST_CHECK_EQUAL(
          5lu,
          cid.value(GeometryID::sensitive_mask, GeometryID::sensitive_shift));
      BOOST_CHECK_EQUAL(
          6lu, cid.value(GeometryID::channel_mask, GeometryID::channel_shift));
    }
  }

}  // end of namespace Test
}  // end of namespace Acts

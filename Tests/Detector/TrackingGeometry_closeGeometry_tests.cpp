///  This file is part of the ACTS project.
///
///  Copyright (C) 2016 ACTS project team
///
///  This Source Code Form is subject to the terms of the Mozilla Public
///  License, v. 2.0. If a copy of the MPL was not distributed with this
///  file, You can obtain one at http:/// mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE GeometryID Tests
#include <boost/test/included/unit_test.hpp>
#include "ACTS/Utilities/Units.hpp"
#include "GeometryCreation.hpp"

namespace Acts {

namespace Test {

  ///  create three cylinder surfaces
  ///  the surface radius (will also be the layer radius)
  double iv_surfaceHalfLengthZ = 50 * Acts::units::_mm;
  double iv_surfaceRadius      = 25. * Acts::units::_mm;
  double iv_surfaceRstagger    = 5. * Acts::units::_mm;
  double iv_surfaceZoverlap    = 10. * Acts::units::_mm;
  double iv_layerEnvelope      = 0.5 * Acts::units::_mm;
  double iv_volumeEnvelope     = 10. * Acts::units::_mm;
  double iv_volumeRadius       = iv_surfaceRadius + 0.5 * iv_surfaceRstagger
      + iv_layerEnvelope + iv_volumeEnvelope;

  ///  the surface radius (will also be the layer radius)
  double ov_surfaceHalfLengthZ = 50. * Acts::units::_mm;
  double ov_surfaceRadius      = 100. * Acts::units::_mm;
  double ov_surfaceRstagger    = 5. * Acts::units::_mm;
  double ov_surfaceZoverlap    = 10. * Acts::units::_mm;
  double ov_layerEnvelope      = 0.5 * Acts::units::_mm;
  double ov_volumeEnvelope     = 10. * Acts::units::_mm;
  double ov_volumeRadius       = ov_surfaceRadius + 0.5 * ov_surfaceRstagger
      + ov_layerEnvelope + ov_volumeEnvelope;

  ///  inner volume
  auto iVolume = constructCylinderVolume(iv_surfaceHalfLengthZ,
                                         iv_surfaceRadius,
                                         iv_surfaceRstagger,
                                         iv_surfaceZoverlap,
                                         iv_layerEnvelope,
                                         iv_volumeEnvelope,
                                         0.,
                                         iv_volumeRadius,
                                         "InnerVolume");
  ///  outer volume
  auto oVolume = constructCylinderVolume(ov_surfaceHalfLengthZ,
                                         ov_surfaceRadius,
                                         ov_surfaceRstagger,
                                         ov_surfaceZoverlap,
                                         ov_layerEnvelope,
                                         ov_volumeEnvelope,
                                         iv_volumeRadius,
                                         ov_volumeRadius,
                                         "OuterVolume");

  // now create the container volume
  double ov_volumeHalfZ
      = (4 * ov_surfaceHalfLengthZ - ov_surfaceZoverlap) + ov_volumeEnvelope;
  auto hVolume = constructContainerVolume(iVolume,
                                          oVolume,
                                          ov_volumeRadius,
                                          ov_volumeHalfZ,
                                          "Container");

  // creating a TrackingGeometry closs the geometry
  TrackingGeometry tGeometry(hVolume);
  // get the world back
  auto world = tGeometry.highestTrackingVolume();

  BOOST_AUTO_TEST_CASE(GeometryID_closeGeometry_test)
  {
    std::cout << "-> Testing Container Volume '" << world->volumeName()
              << std::endl;
    ///  world and containers have to have geoID values 0
    BOOST_CHECK_EQUAL(0, world->geoID().value());
    /// check the boundaries of the world, should also be 0
    std::cout << "-[o] Testing Boundary Surfaces." << std::endl;
    for (auto wbsf : world->boundarySurfaces()) {
      BOOST_CHECK_EQUAL(0, wbsf->surfaceRepresentation().geoID().value());
    }
    /// go through the sub volumes
    geo_id_value vol_id = 0;
    for (auto cVol : hVolume->confinedVolumes()->arrayObjects()) {
      /// let's check everything is set to 0
      std::cout << "--> Testing Volume '" << cVol->volumeName() << std::endl;
      geo_id_value c_vol_id = cVol->geoID().value(GeometryID::volume_mask,
                                                  GeometryID::volume_shift);
      BOOST_CHECK_EQUAL(++vol_id, c_vol_id);
      // check the boundary surfaces
      geo_id_value bsurface_id = 0;
      std::cout << "--[o] Testing Boundary Surfaces." << std::endl;
      for (auto bSf : cVol->boundarySurfaces()) {
        // check the bsurface volume id
        geo_id_value bs_vol_id = bSf->surfaceRepresentation().geoID().value(
            GeometryID::volume_mask, GeometryID::volume_shift);
        geo_id_value bs_bsf_id = bSf->surfaceRepresentation().geoID().value(
            GeometryID::boundary_mask, GeometryID::boundary_shift);
        BOOST_CHECK_EQUAL(++bsurface_id, bs_bsf_id);
        BOOST_CHECK_EQUAL(c_vol_id, bs_vol_id);
      }
      // testing the layer and it's approach surfaces
      geo_id_value layer_id = 0;
      std::cout << "--[o] Testing Layer, the approach and sub surfaces."
                << std::endl;
      for (auto lay : cVol->confinedLayers()->arrayObjects()) {
        // check the layer volume id and layer layer id
        geo_id_value lay_vol_id = lay->geoID().value(GeometryID::volume_mask,
                                                     GeometryID::volume_shift);
        geo_id_value lay_lay_id = lay->geoID().value(GeometryID::layer_mask,
                                                     GeometryID::layer_shift);
        BOOST_CHECK_EQUAL(++layer_id, lay_lay_id);
        BOOST_CHECK_EQUAL(c_vol_id, lay_vol_id);
        // test the layer approach surfaces
        geo_id_value asurface_id = 0;
        for (auto asf : lay->approachDescriptor()->containedSurfaces()) {
          // check the approach volume id, approach layer id, approach approach
          // id
          geo_id_value asf_vol_id = asf->geoID().value(
              GeometryID::volume_mask, GeometryID::volume_shift);
          geo_id_value asf_lay_id = asf->geoID().value(GeometryID::layer_mask,
                                                       GeometryID::layer_shift);
          geo_id_value asf_asf_id = asf->geoID().value(
              GeometryID::approach_mask, GeometryID::approach_shift);
          BOOST_CHECK_EQUAL(layer_id, asf_lay_id);
          BOOST_CHECK_EQUAL(c_vol_id, asf_vol_id);
          BOOST_CHECK_EQUAL(++asurface_id, asf_asf_id);
        }
        // get the sub surfaces
        geo_id_value ssurface_id = 0;
        for (auto ssf : lay->surfaceArray()->arrayObjects()) {
          // check the approach volume id, approach layer id, approach approach
          // id
          geo_id_value ssf_vol_id = ssf->geoID().value(
              GeometryID::volume_mask, GeometryID::volume_shift);
          geo_id_value ssf_lay_id = ssf->geoID().value(GeometryID::layer_mask,
                                                       GeometryID::layer_shift);
          geo_id_value ssf_ssf_id = ssf->geoID().value(
              GeometryID::sensitive_mask, GeometryID::sensitive_shift);
          BOOST_CHECK_EQUAL(layer_id, ssf_lay_id);
          BOOST_CHECK_EQUAL(c_vol_id, ssf_vol_id);
          BOOST_CHECK_EQUAL(++ssurface_id, ssf_ssf_id);
        }
      }
    }
  }
}  //  end of namespace Test
}  //  end of namespace Acts

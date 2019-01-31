// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE GeometryID Tests

#include <boost/test/included/unit_test.hpp>
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"
#include "TrackingVolumeCreation.hpp"

namespace Acts {

namespace Test {

  /// we test a two-level hierarchy
  /// every deeper level hierarchy is a a derivate of this
  ///
  /// WorldVolume   with volumeID       == 1
  /// - InnerVolume with volumeID       == 2
  /// -- InnerInnerVolume with volumeID == 3
  /// -- InnerOuterVolume with volumeID == 4
  /// - OuterVolume with volumeID       == 5

  // sensitive surface definitions
  double surfaceHalfLengthZ = 50 * Acts::units::_mm;
  double surfaceRstagger    = 5. * Acts::units::_mm;
  double surfaceZoverlap    = 10. * Acts::units::_mm;
  double layerEnvelope      = 0.5 * Acts::units::_mm;
  double volumeEnvelope     = 10. * Acts::units::_mm;

  // inner inner volume definitions
  double iiv_surfaceRadius = 25. * Acts::units::_mm;
  double iiv_volumeRadius  = iiv_surfaceRadius + 0.5 * surfaceRstagger
      + layerEnvelope + volumeEnvelope;

  ///  inner outer volume defininitions
  double iov_surfaceRadius = 100. * Acts::units::_mm;
  double iov_volumeRadius  = iov_surfaceRadius + 0.5 * surfaceRstagger
      + layerEnvelope + volumeEnvelope;

  ///  inner inner volume
  auto iiVolume = constructCylinderVolume(surfaceHalfLengthZ,
                                          iiv_surfaceRadius,
                                          surfaceRstagger,
                                          surfaceZoverlap,
                                          layerEnvelope,
                                          volumeEnvelope,
                                          0.,
                                          iiv_volumeRadius,
                                          "InnerInnerVolume");
  ///  inner outer volume
  auto ioVolume = constructCylinderVolume(surfaceHalfLengthZ,
                                          iov_surfaceRadius,
                                          surfaceRstagger,
                                          surfaceZoverlap,
                                          layerEnvelope,
                                          volumeEnvelope,
                                          iiv_volumeRadius,
                                          iov_volumeRadius,
                                          "InnerOuterVolume");

  // now create the Inner Container volume
  double volumeHalfZ
      = (4 * surfaceHalfLengthZ - surfaceZoverlap) + volumeEnvelope;
  /// the inner volume
  auto iVolume = constructContainerVolume(iiVolume,
                                          ioVolume,
                                          iov_volumeRadius,
                                          volumeHalfZ,
                                          "InnerVolume");

  // outer volume definitions
  double ov_surfaceRadius = 150. * Acts::units::_mm;
  double ov_volumeRadius  = ov_surfaceRadius + 0.5 * surfaceRstagger
      + layerEnvelope + volumeEnvelope;

  ///  inner outer volume
  auto oVolume = constructCylinderVolume(surfaceHalfLengthZ,
                                         ov_surfaceRadius,
                                         surfaceRstagger,
                                         surfaceZoverlap,
                                         layerEnvelope,
                                         volumeEnvelope,
                                         iov_volumeRadius,
                                         ov_volumeRadius,
                                         "OuterVolume");
  /// the inner volume
  auto volume = constructContainerVolume(iVolume,
                                         oVolume,
                                         ov_volumeRadius,
                                         volumeHalfZ,
                                         "WorldVolume");

  // creating a TrackingGeometry
  // -> closs the geometry, this should set the GeometryID
  TrackingGeometry tGeometry(volume);
  // get the world back
  auto world = tGeometry.highestTrackingVolume();

  BOOST_AUTO_TEST_CASE(GeometryID_closeGeometry_test)
  {

    // the lambda for checking
    auto check_vol = [](const TrackingVolume& vol, geo_id_value geoid) -> void {
      // check the geometry id of the volume
      BOOST_CHECK_EQUAL(geoid, vol.geoID().value(GeometryID::volume_mask));
      // check the geometry id of all boundary surfaces of the volume
      // - this is strictly only possible when glueing if OFF
      geo_id_value bsurface_id = 0;
      for (auto bSf : vol.boundarySurfaces()) {
        // check the bsurface volume id
        geo_id_value bs_vol_id = bSf->surfaceRepresentation().geoID().value(
            GeometryID::volume_mask);
        BOOST_CHECK_EQUAL(geoid, bs_vol_id);
        // check the bsurface boundary id
        geo_id_value bs_bsf_id = bSf->surfaceRepresentation().geoID().value(
            GeometryID::boundary_mask);
        BOOST_CHECK_EQUAL(++bsurface_id, bs_bsf_id);
      }
      // testing the layer and it's approach surfaces
      if (vol.confinedLayers() != nullptr) {
        // layers start are counted from 1 - n
        geo_id_value layer_id = 0;
        for (auto lay : vol.confinedLayers()->arrayObjects()) {
          // check the layer volume id and layer layer id
          geo_id_value lay_vol_id = lay->geoID().value(GeometryID::volume_mask);
          geo_id_value lay_lay_id = lay->geoID().value(GeometryID::layer_mask);
          BOOST_CHECK_EQUAL(++layer_id, lay_lay_id);
          BOOST_CHECK_EQUAL(geoid, lay_vol_id);
          // test the layer approach surfaces
          if (lay->approachDescriptor() != nullptr) {
            // approach surfacesare counted from 1 - n
            geo_id_value asurface_id = 0;
            for (auto asf : lay->approachDescriptor()->containedSurfaces()) {
              // check the approach volume id, approach layer id, approach
              // approach
              // id
              geo_id_value asf_vol_id
                  = asf->geoID().value(GeometryID::volume_mask);
              geo_id_value asf_lay_id
                  = asf->geoID().value(GeometryID::layer_mask);
              geo_id_value asf_asf_id
                  = asf->geoID().value(GeometryID::approach_mask);
              BOOST_CHECK_EQUAL(layer_id, asf_lay_id);
              BOOST_CHECK_EQUAL(geoid, asf_vol_id);
              BOOST_CHECK_EQUAL(++asurface_id, asf_asf_id);
            }
          }
          // test the sensitive surfaces
          if (lay->surfaceArray() != nullptr) {
            // sensitive surfaces are counted from 1 - n
            geo_id_value ssurface_id = 0;
            for (auto ssf : lay->surfaceArray()->surfaces()) {
              // check the approach volume id, approach layer id, approach
              // approach
              // id
              geo_id_value ssf_vol_id
                  = ssf->geoID().value(GeometryID::volume_mask);
              geo_id_value ssf_lay_id
                  = ssf->geoID().value(GeometryID::layer_mask);
              geo_id_value ssf_ssf_id
                  = ssf->geoID().value(GeometryID::sensitive_mask);
              BOOST_CHECK_EQUAL(layer_id, ssf_lay_id);
              BOOST_CHECK_EQUAL(geoid, ssf_vol_id);
              BOOST_CHECK_EQUAL(++ssurface_id, ssf_ssf_id);
            }
          }
        }
      }
    };

    // get the two volumes the world is built of
    auto ioVolumes = world->confinedVolumes()->arrayObjects();
    // check the size - has to be two volumes
    BOOST_CHECK_EQUAL(2ul, ioVolumes.size());
    // get the innermost volumes
    auto iioVolumes = ioVolumes[0]->confinedVolumes()->arrayObjects();
    // check the size - has to be two volumes
    BOOST_CHECK_EQUAL(2ul, iioVolumes.size());

    // check the world
    check_vol(*world, 1);
    // - check InnerVolume
    check_vol(*ioVolumes[0], 2);
    // -- check the InnerInnerVolume
    check_vol(*iioVolumes[0], 3);
    // -- check the InenerOuterVolume
    check_vol(*iioVolumes[1], 4);
    // - check the OuterVolume
    check_vol(*ioVolumes[1], 5);
  }

  BOOST_AUTO_TEST_CASE(TrackingGeometry_testVisitSurfaces)
  {
    // this will also cover TrackingVolume::visitSurfaces
    // its a pretty bare bones test, and only asserts that the
    // method is called on the expected number of surfaces
    size_t nSurfaces = 0;
    tGeometry.visitSurfaces([&nSurfaces](const auto*) { nSurfaces++; });

    BOOST_CHECK_EQUAL(nSurfaces, 9);
  }

}  //  end of namespace Test
}  //  end of namespace Acts

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
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Utilities/Units.hpp"
#include "TrackingVolumeCreation.hpp"

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

  BOOST_AUTO_TEST_CASE(GeometryID_innervolume_test)
  {
    BOOST_CHECK_EQUAL(0ul, iVolume->geoID().value());
    // check the boundary surfaces
    for (auto bSf : iVolume->boundarySurfaces()) {
      BOOST_CHECK_EQUAL(0ul, bSf->surfaceRepresentation().geoID().value());
      for (auto lay : iVolume->confinedLayers()->arrayObjects()) {
        BOOST_CHECK_EQUAL(0ul, lay->geoID().value());
        // check the approach surfaces
        for (auto asf : lay->approachDescriptor()->containedSurfaces()) {
          BOOST_CHECK_EQUAL(0ul, asf->geoID().value());
        }
        // check the layer surface array
        for (auto ssf : lay->surfaceArray()->surfaces()) {
          BOOST_CHECK_EQUAL(0ul, ssf->geoID().value());
        }
      }
    }
  }

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

  BOOST_AUTO_TEST_CASE(GeometryID_outervolume_test)
  {
    BOOST_CHECK_EQUAL(0ul, oVolume->geoID().value());
    // check the boundary surfaces
    for (auto bSf : iVolume->boundarySurfaces()) {
      BOOST_CHECK_EQUAL(0ul, bSf->surfaceRepresentation().geoID().value());
      for (auto lay : oVolume->confinedLayers()->arrayObjects()) {
        BOOST_CHECK_EQUAL(0ul, lay->geoID().value());
        // check the approach surfaces
        for (auto asf : lay->approachDescriptor()->containedSurfaces()) {
          BOOST_CHECK_EQUAL(0ul, asf->geoID().value());
        }
        // check the layer surface array
        for (auto ssf : lay->surfaceArray()->surfaces()) {
          BOOST_CHECK_EQUAL(0ul, ssf->geoID().value());
        }
      }
    }
  }
  //
  double ov_volumeHalfZ
      = (4 * ov_surfaceHalfLengthZ - ov_surfaceZoverlap) + ov_volumeEnvelope;
  // now create the container volume
  auto hVolume = constructContainerVolume(iVolume,
                                          oVolume,
                                          ov_volumeRadius,
                                          ov_volumeHalfZ,
                                          "Container");

  ///  pre-check on GeometryID
  BOOST_AUTO_TEST_CASE(GeometryID_containervolume_test)
  {
    ///  let's check that the geometry ID values are all 0
    BOOST_CHECK_EQUAL(0ul, hVolume->geoID().value());
    /// check the boundaries of the hVolume, should also be 0
    for (auto hbsf : hVolume->boundarySurfaces()) {
      BOOST_CHECK_EQUAL(0ul, hbsf->surfaceRepresentation().geoID().value());
    }
    for (auto cVol : hVolume->confinedVolumes()->arrayObjects()) {
      /// let's check everything is set to 0
      BOOST_CHECK_EQUAL(0ul, cVol->geoID().value());
      // check the boundary surfaces
      for (auto bSf : cVol->boundarySurfaces()) {
        BOOST_CHECK_EQUAL(0ul, bSf->surfaceRepresentation().geoID().value());
      }
      for (auto lay : cVol->confinedLayers()->arrayObjects()) {
        BOOST_CHECK_EQUAL(0ul, lay->geoID().value());
        // check the approach surfaces
        for (auto asf : lay->approachDescriptor()->containedSurfaces()) {
          BOOST_CHECK_EQUAL(0ul, asf->geoID().value());
        }
        // check the layer surface array
        for (auto ssf : lay->surfaceArray()->surfaces()) {
          BOOST_CHECK_EQUAL(0ul, ssf->geoID().value());
        }
      }
    }
  }

}  //  end of namespace Test
}  //  end of namespace Acts

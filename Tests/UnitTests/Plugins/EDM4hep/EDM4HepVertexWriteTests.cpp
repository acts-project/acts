// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsPlugins/EDM4hep/EDM4hepUtil.hpp"

#include <edm4hep/VertexCollection.h>

using namespace Acts;
using namespace ActsPlugins;

namespace {
auto gctx = GeometryContext::dangerouslyDefaultConstruct();
}  // namespace

BOOST_AUTO_TEST_SUITE(EDM4hepVertexWriteTest)

BOOST_AUTO_TEST_CASE(WriteVertex) {
  Vertex vertex;
  vertex.setPosition(Vector3(1, 2, 3));
  vertex.setTime(42);
  SquareMatrix4 covariance;
  // clang-format off
  covariance << 1, 2, 3, 4,
                 2, 5, 6, 7,
                 3, 6, 8, 9,
                 4, 7, 9, 10;
  // clang-format on
  vertex.setFullCovariance(covariance);
  vertex.setFitQuality(7, 8);

  edm4hep::VertexCollection vertices;
  auto to = vertices.create();

  EDM4hepUtil::writeVertex(vertex, to);

  BOOST_CHECK_EQUAL(to.getPosition()[0], 1);
  BOOST_CHECK_EQUAL(to.getPosition()[1], 2);
  BOOST_CHECK_EQUAL(to.getPosition()[2], 3);

  BOOST_CHECK_EQUAL(to.getCovMatrix().getValue(edm4hep::FourMomCoords::x,
                                               edm4hep::FourMomCoords::x),
                    covariance(eFreePos0, eFreePos0));
  BOOST_CHECK_EQUAL(to.getCovMatrix().getValue(edm4hep::FourMomCoords::x,
                                               edm4hep::FourMomCoords::y),
                    covariance(eFreePos0, eFreePos1));
  BOOST_CHECK_EQUAL(to.getCovMatrix().getValue(edm4hep::FourMomCoords::x,
                                               edm4hep::FourMomCoords::z),
                    covariance(eFreePos0, eFreePos2));
  BOOST_CHECK_EQUAL(to.getCovMatrix().getValue(edm4hep::FourMomCoords::y,
                                               edm4hep::FourMomCoords::x),
                    covariance(eFreePos1, eFreePos0));
  BOOST_CHECK_EQUAL(to.getCovMatrix().getValue(edm4hep::FourMomCoords::y,
                                               edm4hep::FourMomCoords::y),
                    covariance(eFreePos1, eFreePos1));
  BOOST_CHECK_EQUAL(to.getCovMatrix().getValue(edm4hep::FourMomCoords::y,
                                               edm4hep::FourMomCoords::z),
                    covariance(eFreePos1, eFreePos2));
  BOOST_CHECK_EQUAL(to.getCovMatrix().getValue(edm4hep::FourMomCoords::z,
                                               edm4hep::FourMomCoords::x),
                    covariance(eFreePos2, eFreePos0));
  BOOST_CHECK_EQUAL(to.getCovMatrix().getValue(edm4hep::FourMomCoords::z,
                                               edm4hep::FourMomCoords::y),
                    covariance(eFreePos2, eFreePos1));
  BOOST_CHECK_EQUAL(to.getCovMatrix().getValue(edm4hep::FourMomCoords::z,
                                               edm4hep::FourMomCoords::z),
                    covariance(eFreePos2, eFreePos2));

  if constexpr (EDM4hepUtil::detail::edm4hepVertexHasTime<
                    edm4hep::MutableVertex>) {
    BOOST_CHECK_EQUAL(to.getPosition()[3], 42);
    BOOST_CHECK_EQUAL(to.getCovMatrix().getValue(edm4hep::FourMomCoords::t,
                                                 edm4hep::FourMomCoords::t),
                      covariance(eFreeTime, eFreeTime));
    BOOST_CHECK_EQUAL(to.getCovMatrix().getValue(edm4hep::FourMomCoords::x,
                                                 edm4hep::FourMomCoords::t),
                      covariance(eFreePos0, eFreeTime));
    BOOST_CHECK_EQUAL(to.getCovMatrix().getValue(edm4hep::FourMomCoords::y,
                                                 edm4hep::FourMomCoords::t),
                      covariance(eFreePos1, eFreeTime));
    BOOST_CHECK_EQUAL(to.getCovMatrix().getValue(edm4hep::FourMomCoords::z,
                                                 edm4hep::FourMomCoords::t),
                      covariance(eFreePos2, eFreeTime));
  }
}

BOOST_AUTO_TEST_SUITE_END()

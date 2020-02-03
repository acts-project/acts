// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Radial Bounds Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <limits>
#include <fstream>

#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

namespace Test {
BOOST_AUTO_TEST_SUITE(Surfaces)

double minRadius = 7.2;
double maxRadius = 12.0;
double minPhi = 0.74195;
double maxPhi = 1.33970;

Vector2D offset(-2., 2.);

bool writeObj = false;

static void writeAnnulusDiscObj(std::ofstream& stream, double scalor,
                                unsigned int nSegments,
                                const Acts::Transform3D& transform,
                                const std::vector<Acts::Vector2D>& vertices,
                                const Acts::Vector3D& os = Acts::Vector3D(0.,0.,0.)) {
  unsigned int cvc = 0;

  std::vector<Acts::Vector3D> gVertices;
  // Write the vertices
  for (auto& v : vertices) {
    gVertices.push_back(transform * Acts::Vector3D(v.x(), v.y(), 0.));
  }

  // Complete if there are segments defined for a bow
  if (nSegments > 1) {
    /// Fill the vertices in betwen
    auto fillBow =
        [&](const Acts::Vector3D& first,
            const Acts::Vector3D& second) -> std::vector<Acts::Vector3D> {
      // The completed list of vertices
      std::vector<Acts::Vector3D> completed3D;

      double oseg = 1. / nSegments;
      double phif = Acts::VectorHelpers::phi(first);
      double phis = Acts::VectorHelpers::phi(second);
      double phiD = (phis - phif);

      if (std::abs(phiD) > M_PI and phif * phis < 0.) {
        phiD += phiD < 0. ? 2 * M_PI : -2 * M_PI;
      }

      double rf = Acts::VectorHelpers::perp(first);
      double rs = Acts::VectorHelpers::perp(second);
      double zf = first.z();
      double zs = second.z();

      phiD *= oseg;
      double rD = (rs - rf) * oseg;
      double zD = (zs - zf) * oseg;

      for (unsigned int is = 0; is < nSegments + 1; ++is) {
        double r = rf + is * rD;
        double phi = phif + is * phiD;
        double z = zf + is * zD;
        completed3D.push_back(Acts::Vector3D(r * cos(phi), r * sin(phi), z) -
                              os);
      }
      // Reassing the global points
      return completed3D;
    };
    // Fill the bows
    auto completedBow1 = fillBow(gVertices[0], gVertices[1]);
    auto completedBow2 = fillBow(gVertices[2], gVertices[3]);
    // Clear the vertices
    gVertices = completedBow1;
    gVertices.insert(gVertices.end(), completedBow2.begin(),
                     completedBow2.end());
  }

  // Write the vertices
  for (const auto& gv : gVertices) {
    stream << "v " << scalor * gv.x() << " " << scalor * gv.y() << " "
           << scalor * gv.z() << std::endl;
  }
  // Write the faces
  stream << "f";
  for (unsigned int iv = 0; iv < gVertices.size(); ++iv) {
    stream << " " << cvc + iv + 1;
  }
  stream << std::endl;
}

/// Unit tests for AnnulusBounds constrcuctors
BOOST_AUTO_TEST_CASE(AnnulusBoundsConstruction) {
  //
  /// Test construction with radii and default sector
  auto original = AnnulusBounds(minRadius, maxRadius, minPhi, maxPhi, offset);
  BOOST_CHECK_EQUAL(original.type(), SurfaceBounds::Annulus);

  AnnulusBounds copied(original);
  BOOST_CHECK_EQUAL(copied.type(), SurfaceBounds::Annulus);

  if (writeObj) {
    std::ofstream rOut;
    rOut.open("AnnulusDisc.obj");
    unsigned int nSegments = 72;
    double phiStep = 2 * M_PI / nSegments;
    for (unsigned int is = 0; is < nSegments; ++is) {
      double phi = is * phiStep - M_PI;
      double cphi = std::cos(phi);
      double sphi = std::sin(phi);
      rOut << "v " << minRadius * cphi << " " << minRadius * sphi << " 0."
           << std::endl;
      rOut << "v " << maxRadius * cphi << " " << maxRadius * sphi << " 0."
           << std::endl;
    }
    for (unsigned int il = 1; il < 2 * nSegments; ++il) {
      rOut << "l " << il << " " << il + 2 << std::endl;
    }
    rOut.close();

    std::ofstream abOut;
    abOut.open("AnnulusBoundTests.obj");
    writeAnnulusDiscObj(abOut, 1., 12, Transform3D::Identity(),
                        original.vertices(),
                        Vector3D(offset.x(), offset.y(), 0.));
    abOut << std::endl;
    abOut.close();

    std::ofstream linesOut;
    linesOut.open("AnnulusLines.obj");
    auto writeVertex = [&](const Vector2D& v,
                           const Vector2D& o = Vector2D(0., 0.)) -> void {
      linesOut << "v " << v.x() - o.x() << " " << v.y() - o.y() << " 0"
               << std::endl;
    };

    // write extra lines
    auto vertices = original.vertices();
    Vector2D sideA = vertices[0] - vertices[3];
    Vector2D sideB = vertices[1] - vertices[2];

    writeVertex(-offset);
    writeVertex(vertices[3] - 10 * sideA, offset);
    writeVertex(vertices[3] + 11 * sideA, offset);
    writeVertex(vertices[2] - 10 * sideB, offset);
    writeVertex(vertices[2] + 11 * sideB, offset);
    linesOut << "l 2 3" << std::endl;
    linesOut << "l 4 5" << std::endl;
    linesOut.close();
  }
}

/// Unit tests for AnnulusBounds properties
BOOST_AUTO_TEST_CASE(AnnulusBoundsProperties) {
  /// Test construction with radii and default sector
  AnnulusBounds aBounds(minRadius, maxRadius, minPhi, maxPhi, offset);

  /// Test clone
  auto pClonedAnnulusBounds = aBounds.clone();
  BOOST_CHECK_NE(pClonedAnnulusBounds, nullptr);
  delete pClonedAnnulusBounds;
  //
  /// Test type() (redundant; already used in constructor confirmation)
  BOOST_CHECK_EQUAL(aBounds.type(), SurfaceBounds::Annulus);

  /// Test positions inside/outside
  // - start from cartesian (from test drawing)
  Vector2D inSurfaceXY(7., 7.);
  Vector2D outsideXY1(5., 5.);
  Vector2D outsideXY2(10., 3.);
  Vector2D outsideXY3(10., 10.);
  Vector2D outsideXY4(4., 10.);
  std::vector<Vector2D> testPoints = {inSurfaceXY, outsideXY1, outsideXY2,
                                      outsideXY3, outsideXY4};

  auto toStripFrame = [&](const Vector2D& xy) -> Vector2D {
    auto shifted = xy + offset;
    double r = VectorHelpers::perp(shifted);
    double phi = VectorHelpers::phi(shifted);
    return Vector2D(r, phi);
  };

  BOOST_CHECK(aBounds.inside(toStripFrame(inSurfaceXY), BoundaryCheck(true)));
  BOOST_CHECK(!aBounds.inside(toStripFrame(outsideXY1), BoundaryCheck(true)));
  BOOST_CHECK(!aBounds.inside(toStripFrame(outsideXY2), BoundaryCheck(true)));
  BOOST_CHECK(!aBounds.inside(toStripFrame(outsideXY3), BoundaryCheck(true)));
  BOOST_CHECK(!aBounds.inside(toStripFrame(outsideXY4), BoundaryCheck(true)));

  /// Check radial inside
  BOOST_CHECK(!aBounds.insideRadialBounds(0.5));
  BOOST_CHECK(aBounds.insideRadialBounds(9.));
  BOOST_CHECK(!aBounds.insideRadialBounds(18.));

  /// Check binning value
  // double binningValueR();

  /// Check the value in phi
  // CHECK_CLOSE_ABS(0.5*(minPhi+maxPhi), aBounds.binningValuePhi(), 1e-6);

  /// Test rMin
  BOOST_CHECK_EQUAL(aBounds.rMin(), minRadius);
  //
  /// Test rMax
  BOOST_CHECK_EQUAL(aBounds.rMax(), maxRadius);
  /// Test phiMin
  BOOST_CHECK_EQUAL(aBounds.rMin(), minRadius);
  //
  /// Test phiMax
  BOOST_CHECK_EQUAL(aBounds.rMax(), maxRadius);

  if (writeObj) {
    std::ofstream testOut;
    testOut.open("AnnulusTestPoints.obj");
    for (const auto& vc : testPoints) {
      testOut << "v " << vc.x() << " " << vc.y() << " " << 0 << std::endl;
    }
    testOut.close();
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts

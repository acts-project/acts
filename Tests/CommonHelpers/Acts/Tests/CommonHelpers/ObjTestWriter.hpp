// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <fstream>
#include <vector>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ObjHelper.hpp"

namespace Acts {

using IdentifiedPolyderon = std::tuple<std::string, bool, Polyhedron>;

namespace Test {

struct ObjTestWriter {
  /// Helper method to write sector lines for 2D views
  /// @param name The output file name
  /// @param linaA The first line
  /// @param lineB The second line
  static void writeSectorLinesObj(const std::string& name,
                                  const std::pair<Vector3D, Vector3D>& lineA,
                                  const std::pair<Vector3D, Vector3D>& lineB) {
    std::ofstream ostream;
    ostream.open(name + ".obj");
    ObjHelper objH;
    objH.line(lineA.first, lineA.second);
    objH.line(lineB.first, lineB.second);
    objH.write(ostream);
    ostream.close();
  }

  /// Helper method to write sector planes for two dimensional sectors
  /// (symmetric)
  /// @param name The output file name
  /// @param phiSec The opening angle
  /// @param phiAvg The average phi (= center phi position)
  /// @param hX The half length in X of the sector plane
  /// @param hY the half length in Y of the sector plane
  static void writeSectorPlanesObj(const std::string& name, double phiSec,
                                   double phiAvg, double hX, double hY) {
    // Construct the helper planes for sectoral building
    auto sectorBounds = std::make_shared<RectangleBounds>(hX, hY);

    Vector3D helperColX(0., 0., 1.);
    Vector3D helperColY(1., 0., 0.);
    Vector3D helperColZ(0., 1., 0.);
    RotationMatrix3D helperRotation;
    helperRotation.col(0) = helperColX;
    helperRotation.col(1) = helperColY;
    helperRotation.col(2) = helperColZ;
    // curvilinear surfaces are boundless
    Transform3D helperTransform{helperRotation};

    auto sectorTransformM = std::make_shared<Transform3D>(helperTransform);
    sectorTransformM->prerotate(AngleAxis3D(phiAvg - phiSec, helperColX));

    auto sectorTransformP = std::make_shared<Transform3D>(helperTransform);
    sectorTransformP->prerotate(AngleAxis3D(phiAvg + phiSec, helperColX));

    auto sectorPlaneM =
        Surface::makeShared<PlaneSurface>(sectorTransformM, sectorBounds);

    auto sectorPlaneP =
        Surface::makeShared<PlaneSurface>(sectorTransformP, sectorBounds);

    std::ofstream ostream;
    ostream.open(name + ".obj");
    ObjHelper objH;
    sectorPlaneM->polyhedronRepresentation(GeometryContext(), 1).draw(objH);
    sectorPlaneP->polyhedronRepresentation(GeometryContext(), 1).draw(objH);
    objH.write(ostream);
    ostream.close();
  }

  /// Helper method to be called from sub tests
  /// It will draw the polyhedron and create a file, the boolean
  /// steers whether the obj should be triangulated
  /// @param iphs The Identified Polyhedrons (= with name and boolean)
  static void writeObj(const std::vector<IdentifiedPolyderon>& iphs) {
    for (const auto& iph : iphs) {
      std::ofstream ostream;
      ostream.open(std::get<std::string>(iph) + ".obj");
      ObjHelper objH;
      std::get<Polyhedron>(iph).draw(objH, std::get<bool>(iph));
      objH.write(ostream);
      ostream.close();
    }
  }

  /// Helper method to be called from sub tests to write bounding boxes
  /// It will draw the polyhedron and create a file, the boolean
  /// steers whether the obj should be triangulated
  /// @param iphs The Identified Polyhedrons (= with name and boolean)
  static void writeObj(const std::string& name,
                       const Volume::BoundingBox& bBox) {
    std::ofstream ostream;
    ostream.open(name + std::string(".obj"));
    ObjHelper objH;
    bBox.draw(objH);
    objH.write(ostream);
    ostream.close();
  }
};

}  // namespace Test

}  // namespace Acts
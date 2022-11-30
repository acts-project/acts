// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Geant4/Geant4Converters.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4VSolid.hh"

Acts::Transform3 Acts::Geant4AlgebraConverter::transform(
    const G4ThreeVector& trans) {
  Transform3 gTransform = Transform3::Identity();
  Vector3 scaledTrans =
      Vector3(scale * trans[0], scale * trans[1], scale * trans[2]);
  gTransform.pretranslate(scaledTrans);
  return gTransform;
}

Acts::Transform3 Acts::Geant4AlgebraConverter::transform(
    const G4RotationMatrix& rot, const G4ThreeVector& trans) {
  auto gTransform = transform(trans);
  RotationMatrix3 gRot;
  gRot << rot.xx(), rot.xy(), rot.xz(), rot.yx(), rot.yy(), rot.yz(), rot.zx(),
      rot.zy(), rot.zz();
  gTransform.prerotate(gRot);
  return gTransform;
}

std::shared_ptr<Acts::CylinderBounds>
Acts::Geant4ShapeConverter::cylinderBounds(const G4Tubs& g4Tubs) {
  std::array<Acts::ActsScalar, 6u> tArray = {};
  tArray[0u] = static_cast<ActsScalar>(g4Tubs.GetInnerRadius() +
                                       g4Tubs.GetOuterRadius()) *
               0.5;
  tArray[1u] = static_cast<ActsScalar>(g4Tubs.GetZHalfLength());
  tArray[2u] = 0.5 * static_cast<ActsScalar>(g4Tubs.GetDeltaPhiAngle());
  tArray[3u] = static_cast<ActsScalar>(g4Tubs.GetStartPhiAngle());
  return std::make_shared<CylinderBounds>(tArray);
}

std::shared_ptr<Acts::RadialBounds> Acts::Geant4ShapeConverter::radialBounds(
    const G4Tubs& g4Tubs) {
  std::array<ActsScalar, 4u> tArray = {};
  tArray[0u] = static_cast<ActsScalar>(g4Tubs.GetInnerRadius());
  tArray[1u] = static_cast<ActsScalar>(g4Tubs.GetOuterRadius());
  tArray[2u] = 0.5 * static_cast<ActsScalar>(g4Tubs.GetDeltaPhiAngle());
  tArray[3u] = static_cast<ActsScalar>(g4Tubs.GetStartPhiAngle());
  return std::make_shared<RadialBounds>(tArray);
}

std::tuple<std::shared_ptr<Acts::RectangleBounds>, std::array<int, 2u>>
Acts::Geant4ShapeConverter::rectangleBounds(const G4Box& g4Box) {
  std::vector<ActsScalar> hG4XYZ = {
      static_cast<ActsScalar>(g4Box.GetXHalfLength()),
      static_cast<ActsScalar>(g4Box.GetYHalfLength()),
      static_cast<ActsScalar>(g4Box.GetZHalfLength())};

  auto minAt = std::min_element(hG4XYZ.begin(), hG4XYZ.end());
  std::size_t minPos = std::distance(hG4XYZ.begin(), minAt);

  std::array<int, 2u> rAxes = {};
  switch (minPos) {
    case 0: {
      rAxes = {1, 2};
    } break;
    case 1: {
      if (keepAxisOrder) {
        rAxes = {0, -2};  // flip for right-handed
      } else {
        rAxes = {2, 0};  // cylcic positive
      }
    } break;
    case 2: {
      rAxes = {0, 1};
    } break;
  }
  auto rBounds = std::make_shared<RectangleBounds>(hG4XYZ[std::abs(rAxes[0u])],
                                                   hG4XYZ[std::abs(rAxes[1u])]);
  return std::tie(rBounds, rAxes);
}

std::tuple<std::shared_ptr<Acts::TrapezoidBounds>, std::array<int, 2u>>
Acts::Geant4ShapeConverter::trapezoidBounds(const G4Trd& g4Trd) {
  // primary parameters
  ActsScalar hlX0 = static_cast<ActsScalar>(g4Trd.GetXHalfLength1());
  ActsScalar hlX1 = static_cast<ActsScalar>(g4Trd.GetXHalfLength2());
  ActsScalar hlY0 = static_cast<ActsScalar>(g4Trd.GetYHalfLength1());
  ActsScalar hlY1 = static_cast<ActsScalar>(g4Trd.GetYHalfLength2());
  ActsScalar hlZ = static_cast<ActsScalar>(g4Trd.GetZHalfLength());

  std::vector<ActsScalar> dXYZ = {(hlX0 + hlX1) * 0.5, (hlY0 + hlY1) * 0.5,
                                  hlZ};

  auto minAt = std::min_element(dXYZ.begin(), dXYZ.end());
  std::size_t minPos = std::distance(dXYZ.begin(), minAt);

  ActsScalar halfLengthXminY = 0.;
  ActsScalar halfLengthXmaxY = 0.;
  ActsScalar halfLengthY = 0.;

  std::array<int, 2u> rAxes = {};
  switch (minPos) {
    case 0: {
      halfLengthXminY = hlY0;
      halfLengthXmaxY = hlY1;
      halfLengthY = hlZ;
      rAxes = {1, 2};
    } break;
    case 1: {
      halfLengthXminY = hlX0;
      halfLengthXmaxY = hlX1;
      halfLengthY = hlZ;
      rAxes = {0, -2};
    } break;
    case 2: {
      if (std::abs(hlY0 - hlY1) < std::abs(hlX0 - hlX1)) {
        halfLengthXminY = hlX0;
        halfLengthXmaxY = hlX1;
        halfLengthY = (hlY0 + hlY1) * 0.5;
        rAxes = {0, 1};
      } else {
        halfLengthXminY = hlY0;
        halfLengthXmaxY = hlY1;
        halfLengthY = (hlX0 + hlX1) * 0.5;
        rAxes = {-1, 0};
      }
    } break;
  }

  auto tBounds = std::make_shared<TrapezoidBounds>(
      halfLengthXminY, halfLengthXmaxY, halfLengthY);
  return std::tie(tBounds, rAxes);
}

std::tuple<std::shared_ptr<Acts::PlanarBounds>, std::array<int, 2u>>
Acts::Geant4ShapeConverter::planarBounds(const G4VSolid& g4Solid) {

  const G4Box* box = dynamic_cast<const G4Box*>(&g4Solid);
  if (box != nullptr) {
    auto [rBounds, axes] = rectangleBounds(*box);
    return std::tie(rBounds, axes);
  }

  const G4Trd* trd = dynamic_cast<const G4Trd*>(&g4Solid);
  if (trd != nullptr) {
    auto [tBounds, axes] = trapezoidBounds(*trd);
    return std::tie(tBounds, axes);
  }

  std::shared_ptr<Acts::PlanarBounds> pBounds = nullptr;
  std::array<int, 2u> rAxes = {};
  return std::tie(pBounds, rAxes);
}

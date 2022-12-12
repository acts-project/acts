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
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"

Acts::Transform3 Acts::Geant4AlgebraConverter::transform(
    const G4ThreeVector& g4Trans) {
  Transform3 gTransform = Transform3::Identity();
  Vector3 scaledTrans =
      Vector3(scale * g4Trans[0], scale * g4Trans[1], scale * g4Trans[2]);
  gTransform.pretranslate(scaledTrans);
  return gTransform;
}

Acts::Transform3 Acts::Geant4AlgebraConverter::transform(
    const G4RotationMatrix& g4Rot, const G4ThreeVector& g4Trans) {
  // Create the translation
  Vector3 translation(scale * g4Trans[0], scale * g4Trans[1],
                      scale * g4Trans[2]);
  // And the rotation to it
  RotationMatrix3 rotation;
  rotation << g4Rot.xx(), g4Rot.yx(), g4Rot.zx(), g4Rot.xy(), g4Rot.yy(),
      g4Rot.zy(), g4Rot.xz(), g4Rot.yz(), g4Rot.zz();
  Transform3 transform = Transform3::Identity();
  transform.matrix().block(0, 0, 3, 1) = rotation.col(0);
  transform.matrix().block(0, 1, 3, 1) = rotation.col(1);
  transform.matrix().block(0, 2, 3, 1) = rotation.col(2);
  transform.matrix().block(0, 3, 3, 1) = translation;
  return transform;
}

Acts::Transform3 Acts::Geant4AlgebraConverter::transform(
    const G4Transform3D& g4Trf) {
  // Create the translation
  Vector3 translation(scale * g4Trf.dx(), scale * g4Trf.dy(),
                      scale * g4Trf.dz());
  // And the rotation to it
  RotationMatrix3 rotation;
  rotation << g4Trf.xx(), g4Trf.yx(), g4Trf.zx(), g4Trf.xy(), g4Trf.yy(),
      g4Trf.zy(), g4Trf.xz(), g4Trf.yz(), g4Trf.zz();
  Transform3 transform = Transform3::Identity();
  transform.matrix().block(0, 0, 3, 1) = rotation.col(0);
  transform.matrix().block(0, 1, 3, 1) = rotation.col(1);
  transform.matrix().block(0, 2, 3, 1) = rotation.col(2);
  transform.matrix().block(0, 3, 3, 1) = translation;
  return transform;
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

namespace {
Acts::Transform3 axesOriented(const Acts::Transform3& toGlobalOriginal,
                              const std::array<int, 2u>& axes) {
  auto originalRotation = toGlobalOriginal.rotation();
  auto colX = originalRotation.col(std::abs(axes[0u]));
  auto colY = originalRotation.col(std::abs(axes[1u]));
  colX *= std::copysign(1, axes[0u]);
  colY *= std::copysign(1, axes[1u]);
  Acts::Vector3 colZ = colX.cross(colY);

  Acts::Transform3 orientedTransform = Acts::Transform3::Identity();
  orientedTransform.matrix().block(0, 0, 3, 1) = colX;
  orientedTransform.matrix().block(0, 1, 3, 1) = colY;
  orientedTransform.matrix().block(0, 2, 3, 1) = colZ;
  orientedTransform.matrix().block(0, 3, 3, 1) = toGlobalOriginal.translation();

  return orientedTransform;
}
}  // namespace

std::shared_ptr<Acts::Surface> Acts::Geant4PhysicalVolumeConverter::surface(
    const G4VPhysicalVolume& g4PhysVol, const Transform3& toGlobal) {
  // Get the logical volume
  auto g4LogVol = g4PhysVol.GetLogicalVolume();
  auto g4Solid = g4LogVol->GetSolid();

  // Dynamic cast chain & conversion
  auto g4Box = dynamic_cast<const G4Box*>(g4Solid);
  if (g4Box != nullptr) {
    auto [bounds, axes] = Geant4ShapeConverter{}.rectangleBounds(*g4Box);
    auto orientedToGlobal = axesOriented(toGlobal, axes);
    return Acts::Surface::makeShared<PlaneSurface>(orientedToGlobal,
                                                   std::move(bounds));
  }

  auto g4Trd = dynamic_cast<const G4Trd*>(g4Solid);
  if (g4Trd != nullptr) {
    auto [bounds, axes] = Geant4ShapeConverter{}.trapezoidBounds(*g4Trd);
    auto orientedToGlobal = axesOriented(toGlobal, axes);
    return Acts::Surface::makeShared<PlaneSurface>(orientedToGlobal,
                                                   std::move(bounds));
  }

  return nullptr;
}

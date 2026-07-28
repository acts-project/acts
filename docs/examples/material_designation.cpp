// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <memory>

using namespace Acts::UnitLiterals;

/// Designate material on volume faces that is to be filled in by the material
/// mapping. This is the Gen3 equivalent of setting "mapMaterial": true plus a
/// bin count in the Gen1 JSON geometry map.
void exampleDesignateProtoMaterial(Acts::BlueprintNode& parent) {
  //! [Designate Proto Material]
  parent.addMaterial("PixelMaterial", [&](auto& mat) {
    using enum Acts::AxisDirection;
    using enum Acts::AxisBoundaryType;
    using enum Acts::CylinderVolumeBounds::Face;

    // Two proto axes per face: the mapping fills these bins in later.
    mat.configureFace(OuterCylinder, {AxisRPhi, Bound, 20}, {AxisZ, Bound, 20});
    mat.configureFace(NegativeDisc, {AxisR, Bound, 15}, {AxisPhi, Bound, 25});
    mat.configureFace(PositiveDisc, {AxisR, Bound, 15}, {AxisPhi, Bound, 25});

    // The material node wraps exactly one child: the volume it applies to.
    mat.addCylinderContainer("Pixel", Acts::AxisDirection::AxisZ,
                             [&](auto& /*pixel*/) {
                               // ... build the pixel detector here
                             });
  });
  //! [Designate Proto Material]
}

/// Assign known material directly at construction. Nothing is mapped onto these
/// faces: they already carry their final material.
void exampleDesignateHomogeneousMaterial(Acts::BlueprintNode& parent) {
  //! [Designate Homogeneous Material]
  // Beryllium, 0.8 mm thick.
  auto beampipeMaterial =
      std::make_shared<const Acts::HomogeneousSurfaceMaterial>(
          Acts::MaterialSlab(Acts::Material::fromMassDensity(
                                 352.8_mm, 407_mm, 9.012, 4, 1.848_g / 1_cm3),
                             0.8_mm));

  parent.addMaterial("BeampipeMaterial", [&](auto& bpMat) {
    using enum Acts::CylinderVolumeBounds::Face;

    bpMat.configureFace(OuterCylinder, beampipeMaterial);
    bpMat.addStaticVolume(
        Acts::Transform3::Identity(),
        std::make_shared<Acts::CylinderVolumeBounds>(0_mm, 25_mm, 1000_mm),
        "Beampipe");
  });
  //! [Designate Homogeneous Material]
}

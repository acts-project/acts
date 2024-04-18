// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/DetectorMaterialHelper.hpp"

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"

void Acts::Experimental::DetectorMaterialHelper::assignMaterial(
    Detector& detector, const DetectorMaterialMaps& materialMaps) {
  SurfaceMaterialAssigner surfaceMaterialAssigner{materialMaps.first};
  detector.visitMutableSurfaces(surfaceMaterialAssigner);

  VolumeMaterialAssigner volumeMaterialAssigner{materialMaps.second};
  detector.visitMutableVolumes(volumeMaterialAssigner);
}

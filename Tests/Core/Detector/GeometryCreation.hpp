// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Layers/CylinderLayer.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/CylinderVolumeBounds.hpp"

namespace Acts {

///  helper function to create a cylinder
TrackingVolumePtr
constructCylinderVolume(double             surfaceHalfLengthZ,
                        double             surfaceRadius,
                        double             surfaceRstagger,
                        double             surfaceZoverlap,
                        double             layerEnvelope,
                        double             volumeEnvelope,
                        double             innerVolumeR,
                        double             outerVolumeR,
                        const std::string& name)
{
  ///  the surface transforms
  auto sfnPosition
      = Vector3D(0., 0., -3 * surfaceHalfLengthZ - surfaceZoverlap);
  auto sfnTransform
      = std::make_shared<const Transform3D>(Translation3D(sfnPosition));
  auto sfcTransform = nullptr;
  auto sfpPosition = Vector3D(0., 0., 3 * surfaceHalfLengthZ - surfaceZoverlap);
  auto sfpTransform
      = std::make_shared<const Transform3D>(Translation3D(sfpPosition));
  ///  the surfaces
  auto sfn = Surface::makeShared<CylinderSurface>(
      sfnTransform, surfaceRadius - 0.5 * surfaceRstagger, surfaceHalfLengthZ);
  auto sfc = Surface::makeShared<CylinderSurface>(
      sfcTransform, surfaceRadius + 0.5 * surfaceRstagger, surfaceHalfLengthZ);
  auto sfp = Surface::makeShared<CylinderSurface>(
      sfpTransform, surfaceRadius - 0.5 * surfaceRstagger, surfaceHalfLengthZ);

  ///  prepare the surfaces

  ///  make the binned array
  double bUmin = sfnPosition.z() - surfaceHalfLengthZ;
  double bUmax = sfpPosition.z() + surfaceHalfLengthZ;

  std::vector<std::shared_ptr<const Surface>> surfaces_only = {{sfn, sfc, sfp}};
  std::vector<const Surface*>                 surfaces_only_raw
      = {{sfn.get(), sfc.get(), sfp.get()}};

  detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Bound>
       axis(bUmin, bUmax, surfaces_only.size());
  auto g2l = [](const Vector3D& glob) {
    return std::array<double, 1>({{glob.z()}});
  };
  auto l2g
      = [](const std::array<double, 1>& loc) { return Vector3D(0, 0, loc[0]); };
  auto sl = std::make_unique<SurfaceArray::SurfaceGridLookup<decltype(axis)>>(
      g2l, l2g, std::make_tuple(axis));
  sl->fill(surfaces_only_raw);
  auto bArray = std::make_unique<SurfaceArray>(std::move(sl), surfaces_only);

  ///  now create the Layer
  auto layer0bounds
      = std::make_shared<const CylinderBounds>(surfaceRadius, bUmax);
  auto layer0 = CylinderLayer::create(nullptr,
                                      layer0bounds,
                                      std::move(bArray),
                                      surfaceRstagger + 2 * layerEnvelope);
  std::unique_ptr<const LayerArray> layerArray
      = std::make_unique<const BinnedArrayXD<LayerPtr>>(layer0);

  ///  create the volume
  auto volumeBounds = std::make_shared<const CylinderVolumeBounds>(
      innerVolumeR, outerVolumeR, bUmax + volumeEnvelope);
  TrackingVolumePtr volume = TrackingVolume::create(
      nullptr, volumeBounds, nullptr, std::move(layerArray), {}, {}, {}, name);
  ///  return the volume
  return volume;
}

///  helper function to create a container
MutableTrackingVolumePtr
constructContainerVolume(TrackingVolumePtr  iVolume,
                         TrackingVolumePtr  oVolume,
                         double             hVolumeRadius,
                         double             hVolumeHalflength,
                         const std::string& name)
{
  ///  create the volume array
  using VAP                = std::pair<TrackingVolumePtr, Vector3D>;
  std::vector<VAP> volumes = {{iVolume, iVolume->binningPosition(binR)},
                              {oVolume, oVolume->binningPosition(binR)}};
  ///  the bounds for the container
  auto hVolumeBounds = std::make_shared<const CylinderVolumeBounds>(
      0., hVolumeRadius, hVolumeHalflength);
  ///  create the BinUtility & the BinnedArray
  auto vUtility = std::make_unique<const BinUtility>(
      volumes.size(), 0., hVolumeRadius, open, binR);
  std::shared_ptr<const TrackingVolumeArray> vArray
      = std::make_shared<const BinnedArrayXD<TrackingVolumePtr>>(
          volumes, std::move(vUtility));
  ///  create the container volume
  auto hVolume = TrackingVolume::create(nullptr, hVolumeBounds, vArray, name);
  // return the container
  return hVolume;
}
}

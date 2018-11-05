// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GenericApproachDescriptor.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <algorithm>
#include "Acts/Utilities/ApproachDescriptor.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Volumes/BoundarySurfaceT.hpp"

namespace Acts {

/// @class GenericApproachDescriptor
///
/// Class to decide and return which approaching surface to be taken,
/// it's a generic descriptor for n surfaces
///
/// It is templated in order to allow for BoundarySurfaces from
/// representing volumes of layers to be re-used

class GenericApproachDescriptor : public ApproachDescriptor
{
public:
  /// A generic approach descriptor for new Acts::Surface objects
  /// passing ownership
  ///
  /// @param aSurfaces are the approach surfaces
  GenericApproachDescriptor(
      std::vector<std::shared_ptr<const Surface>> aSurfaces)
    : ApproachDescriptor(), m_surfaces(std::move(aSurfaces)), m_surfaceCache()
  {
    m_surfaceCache = unpack_shared_vector(m_surfaces);
  }

  /// A generic approach descriptor with n surfaces to test
  ~GenericApproachDescriptor() override = default;

  /// register the Layer to the surfaces
  ///
  /// @param lay is the layer to be registerd
  void
  registerLayer(const Layer& lay) override;

  /// get the compatible surfaces
  ///
  /// @param gpos is the global position to start the approach from
  /// @param gdir is the momentum vector
  /// @param bcheck is the boundary check prescription
  /// @param corrfnc is an noption correction function
  ///
  /// @return : a boolean indicating if an actual intersection had been tried
  ObjectIntersection<Surface>
  approachSurface(const Vector3D&      gpos,
                  const Vector3D&      gdir,
                  NavigationDirection  navDir,
                  const BoundaryCheck& bcheck,
                  CorrFnc              corrfnc = nullptr) const override;

  /// return all contained surfaces of this approach descriptor
  const std::vector<const Surface*>&
  containedSurfaces() const override;

  /// Non-const version
  std::vector<const Surface*>&
  containedSurfaces() override;

private:
  /// approach surfaces with ownership control
  std::vector<std::shared_ptr<const Surface>> m_surfaces;

  /// the surface container cache
  ///
  /// We will need to mutate those surfaces in registerLayer, but the C++ type
  /// system has no const-correct way of expressing this constraint.
  ///
  std::vector<const Surface*> m_surfaceCache;
};

}  // namespace Acts
// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ConeLayer.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <algorithm>
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class ConeBounds;
class ApproachDescriptor;

/// @class ConeLayer
///
/// Class to describe a conical detector layer for tracking, it inhertis from
/// both, Layer base class and ConeSurface class
class ConeLayer : virtual public ConeSurface, public Layer {
 public:
  /// Factory for shared layer
  ///
  /// @param transform is the 3D transform that poisitions the layer in 3D frame
  /// @param cbounds is the conical bound description
  /// @param surfaceArray is the array of sensitive surfaces
  /// @param thickness is the layer thickness along the normal axis
  /// @param ad is the approach descriptor for navigation towards the layer
  /// @param laytyp is the layer type
  ///
  /// @todo chage od and ad to unique_ptr
  ///
  /// @return is a shared pointer to a layer
  static MutableLayerPtr create(
      std::shared_ptr<const Transform3D> transform,
      std::shared_ptr<const ConeBounds> cbounds,
      std::unique_ptr<SurfaceArray> surfaceArray, double thickness = 0.,
      std::unique_ptr<ApproachDescriptor> ad = nullptr,
      LayerType laytyp = Acts::active) {
    return MutableLayerPtr(new ConeLayer(
        std::move(transform), std::move(cbounds), std::move(surfaceArray),
        thickness, std::move(ad), laytyp));
  }

  /// Default Constructor - delete
  ConeLayer() = delete;

  /// Copy constructor of ConeLayer - delete
  ConeLayer(const ConeLayer& cla) = delete;

  /// Assignment operator for ConeLayers - delete
  ConeLayer& operator=(const ConeLayer&) = delete;

  /// Destructor
  ~ConeLayer() override = default;

  /// Transforms the layer into a Surface representation for extrapolation
  const ConeSurface& surfaceRepresentation() const override;

  // Non-const version
  ConeSurface& surfaceRepresentation() override;

 protected:
  /// Private constructor with arguments
  ///
  /// @param transform is the 3D transform that poisitions the layer in 3D frame
  /// @param cbounds is the conical bound description
  /// @param surfaceArray is the array of sensitive surfaces
  /// @param thickness is the layer thickness along the normal axis
  /// @param ade is the approach descriptor for navigation towards the layer
  /// @param laytyp is the layer type
  ///
  /// @todo chage od and ad to unique_ptr
  ConeLayer(std::shared_ptr<const Transform3D> transform,
            std::shared_ptr<const ConeBounds> cbounds,
            std::unique_ptr<SurfaceArray> surfaceArray, double thickness = 0.,
            std::unique_ptr<ApproachDescriptor> ade = nullptr,
            LayerType laytyp = Acts::active);

  /// Private copy constructor with shift, called by create(args*)
  ///
  /// @param cla is the source clone layer
  /// @param shift is the additional shift applied after copying
  ConeLayer(const ConeLayer& cla, const Transform3D& shift);
};

}  // namespace Acts
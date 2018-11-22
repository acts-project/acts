// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiscLayer.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <algorithm>
#include "Acts/Layers/Layer.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"

namespace Acts {

class DiscBounds;
class ApproachDescriptor;

/// @class DiscLayer
///
/// Class to describe a disc-like detector layer for tracking,
/// it inhertis from both, Layer base class
/// and DiscSurface class

class DiscLayer : virtual public DiscSurface, public Layer
{
public:
  ///  Factory constructor with DiscSurface components
  ///
  /// @param transform is the transform to place the layer in the 3D frame
  /// @param dbounds are the disc bounds that describe the layer dimensions
  /// @param surfaceArray is the array of sensitive surfaces
  /// @param thickness is the layer thickness (along the normal vector)
  /// @param ad is the approach descriptor that provides the approach surface
  /// @param laytyp is the layer type
  ///
  /// @todo move ApproachDescriptor to unqique_ptr
  ///
  /// @return a sharted pointer to the new layer
  static MutableLayerPtr
  create(const std::shared_ptr<const Transform3D>& transform,
         const std::shared_ptr<const DiscBounds>&  dbounds,
         std::unique_ptr<SurfaceArray>             surfaceArray = nullptr,
         double                                    thickness    = 0.,
         std::unique_ptr<ApproachDescriptor>       ad           = nullptr,
         LayerType                                 laytyp       = Acts::passive)
  {
    return MutableLayerPtr(new DiscLayer(transform,
                                         dbounds,
                                         std::move(surfaceArray),
                                         thickness,
                                         std::move(ad),
                                         laytyp));
  }

  /// Factory for shared Layer pointer, that accepts @c variant_data
  /// @param vardata The data to build from
  static MutableLayerPtr
  create(const variant_data& vardata);

  /// Default Constructor
  DiscLayer() = delete;

  /// Copy constructor of DiscLayer - deleted
  DiscLayer(const DiscLayer& cla) = delete;

  /// Assignment operator for DiscLayers - deleted
  DiscLayer&
  operator=(const DiscLayer&)
      = delete;

  /// Destructor
  ~DiscLayer() override = default;

  /// Transforms the layer into a Surface representation for extrapolation
  /// @return This method returns a surface reference
  const DiscSurface&
  surfaceRepresentation() const override;

  // Non-const version
  DiscSurface&
  surfaceRepresentation() override;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  variant_data
  toVariantData() const override;

private:
  /// build approach surfaces
  void
  buildApproachDescriptor();

protected:
  // Constructor with DiscSurface components and pointer to SurfaceArray
  ///
  /// @param transform is the transform to place the layer in the 3D frame
  /// @param dbounds are the disc bounds that describe the layer dimensions
  /// @param surfaceArray is the array of sensitive surfaces
  /// @param thickness is the layer thickness (along the normal vector)
  /// @param ad is the approach descriptor that provides the approach surface
  /// @param laytyp is the layer taype
  DiscLayer(const std::shared_ptr<const Transform3D>& transform,
            const std::shared_ptr<const DiscBounds>&  dbounds,
            std::unique_ptr<SurfaceArray>             surfaceArray = nullptr,
            double                                    thickness    = 0.,
            std::unique_ptr<ApproachDescriptor>       ades         = nullptr,
            LayerType                                 laytyp = Acts::active);

  /// Copy constructor with shift
  DiscLayer(const DiscLayer& cla, const Transform3D& tr);
};

}  // namespace
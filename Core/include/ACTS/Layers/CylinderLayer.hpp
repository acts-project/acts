// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderLayer.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_LAYERS_CYLINDERLAYER_H
#define ACTS_LAYERS_CYLINDERLAYER_H

#include <algorithm>
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/InstanceFactory.hpp"
#include "ACTS/Utilities/ThrowAssert.hpp"
#include "ACTS/Utilities/VariantData.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"

namespace Acts {

class CylinderBounds;
class SurfaceMaterial;
class ApproachDescriptor;

/// @class CylinderLayer
///
/// Class to describe a cylindrical detector layer for tracking, it inhertis
/// from
/// both,
/// Layer base class and CylinderSurface class
///
class CylinderLayer : public CylinderSurface, public Layer
{
public:
  /// Factory for shared Layer pointer
  /// create a shared, fully deployed CylinderLayer
  ///
  /// @param transform is the 3D transform that places the layer in 3D space
  /// @param cbounds are the cylindrical bounds of the layer
  /// @param surfaceArray is the Binned Array that holds the sensitive surfaces
  /// @param thickness is the layer thickness (along the normal)
  /// @param ad is the approach descriptor for approaching the layer
  /// @param laytyp is the layer type
  ///
  /// @todo ApproachDescriptor to unique_ptr
  ///
  /// @return The return object is a shared poiter to the layer.
  static MutableLayerPtr
  create(std::shared_ptr<const Transform3D>    transform,
         std::shared_ptr<const CylinderBounds> cbounds,
         std::unique_ptr<SurfaceArray>         surfaceArray = nullptr,
         double                                thickness    = 0.,
         std::unique_ptr<ApproachDescriptor>   ad           = nullptr,
         LayerType                             laytyp       = passive)
  {
    return MutableLayerPtr(new CylinderLayer(transform,
                                             cbounds,
                                             std::move(surfaceArray),
                                             thickness,
                                             std::move(ad),
                                             laytyp));
  }

  static MutableLayerPtr
  create(const variant_data& data_)
  {
    throw_assert(data_.which() == 4, "Variant data must be map");
    const variant_map& data = boost::get<variant_map>(data_);
    std::string        type = data.get<std::string>("type");
    throw_assert(type == "CylinderLayer", "Type must be CylinderLayer");

    variant_map payload = data.get<variant_map>("payload");

    auto trf = std::make_shared<const Transform3D>(
        from_variant<Transform3D>(payload.get<variant_map>("transform")));

    LayerType laytyp = passive;

    double thickness = payload.get<double>("thickness");

    InstanceFactory    factory;
    const variant_map& var_bounds = payload.get<variant_map>("cylinder_bounds");

    auto cbounds = std::dynamic_pointer_cast<const CylinderBounds>(
        factory.surfaceBounds(var_bounds.get<std::string>("type"), var_bounds));

    return MutableLayerPtr(new CylinderLayer(trf,
                                             cbounds,
                                             nullptr,  // SurfaceArray
                                             thickness,
                                             nullptr,  // std::move(ad),
                                             laytyp));
  }

  /// Copy constructor - deleted
  CylinderLayer(const CylinderLayer& cla) = delete;

  /// Assignment operator for CylinderLayers - deleted
  CylinderLayer&
  operator=(const CylinderLayer&)
      = delete;

  /// Default Constructor
  CylinderLayer() = delete;

  /// Destructor
  virtual ~CylinderLayer() {}

  /// Transforms the layer into a Surface representation
  /// This is for positioning and extrapolation
  const CylinderSurface&
  surfaceRepresentation() const override;

  // Non-const version
  CylinderSurface&
  surfaceRepresentation() override;

  variant_data
  toVariantData() const
  {
    using namespace std::string_literals;
    variant_map payload;

    payload["transform"] = to_variant(*m_transform);

    // we need to recover the bounds
    const AbstractVolume* absVol = representingVolume();
    auto                  cvBounds
        = dynamic_cast<const CylinderVolumeBounds*>(&absVol->volumeBounds());

    double cylR = cvBounds->innerRadius() + 0.5 * thickness();
    double hlZ  = cvBounds->halflengthZ();

    CylinderBounds cylBounds(cylR, hlZ);

    const variant_data bounds  = cylBounds.toVariantData();
    payload["cylinder_bounds"] = bounds;

    payload["thickness"] = thickness();

    // payload["surface_array"] = m_surfaceArray->toVariantData();
    payload["surface_array"] = 0;

    variant_map data;
    data["type"]    = "CylinderLayer"s;
    data["payload"] = payload;

    return data;
  }

private:
  /// build approach surfaces */
  void
  buildApproachDescriptor();

protected:
  /// Private constructor for CylinderLayer, called by create(args*) factory
  ///
  /// @param transform is the 3D transform that places the layer in 3D space
  /// @param cbounds are the cylindrical bounds of the layer
  /// @param surfaceArray is the Binned Array that holds the sensitive surfaces
  /// @param thickness is the layer thickness (along the normal)
  /// @param ad is the approach descriptor for approaching the layer
  /// @param laytyp is the layer type
  /// @todo change ApproachDescriptor to unique_ptr
  ///
  /// @return The return object is a shared poiter to the layer.
  CylinderLayer(std::shared_ptr<const Transform3D>    transform,
                std::shared_ptr<const CylinderBounds> cbounds,
                std::unique_ptr<SurfaceArray>         surfaceArray = nullptr,
                double                                thickness    = 0.,
                std::unique_ptr<ApproachDescriptor>   ad           = nullptr,
                LayerType                             laytyp       = passive);

  /// Private copy constructor with shift, called by create(args*)
  ///
  /// @param cla is the source cylinder layer for the copy
  /// @param shift is the additional transform applied after cloning
  ///
  /// @return The return object is a shared pointer to the layer.
  CylinderLayer(const CylinderLayer& cla, const Transform3D& shift);
};

}  // namespace Acts

#endif  // ACTS_LAYERS_CYLINDERLAYER_H

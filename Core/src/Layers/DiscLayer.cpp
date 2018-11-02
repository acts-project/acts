// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiscLayer.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Layers/DiscLayer.hpp"
#include "Acts/Layers/GenericApproachDescriptor.hpp"
#include "Acts/Layers/Layer.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantData.hpp"
#include "Acts/Volumes/AbstractVolume.hpp"
#include "Acts/Volumes/BoundarySurfaceFace.hpp"
#include "Acts/Volumes/CylinderVolumeBounds.hpp"

using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::perp;

Acts::DiscLayer::DiscLayer(const std::shared_ptr<const Transform3D>& transform,
                           const std::shared_ptr<const DiscBounds>&  dbounds,
                           std::unique_ptr<SurfaceArray>       surfaceArray,
                           double                              thickness,
                           std::unique_ptr<ApproachDescriptor> ades,
                           LayerType                           laytyp)
  : DiscSurface(transform, dbounds)
  , Layer(std::move(surfaceArray), thickness, std::move(ades), laytyp)
{
  // create the representing volume
  const RadialBounds* rBounds
      = dynamic_cast<const RadialBounds*>(dbounds.get());
  if (rBounds != nullptr) {
    // @todo make a trapezoidal volume when you have DiscTrapezoidalBounds
    CylinderVolumeBounds* cvBounds = new CylinderVolumeBounds(
        rBounds->rMin(), rBounds->rMax(), 0.5 * thickness);
    m_representingVolume
        = new AbstractVolume(transform, VolumeBoundsPtr(cvBounds));
  }
  // associate the layer to the layer surface itself
  DiscSurface::associateLayer(*this);
  // build an approach descriptor if none provided
  if (!m_approachDescriptor && m_surfaceArray) {
    buildApproachDescriptor();
  }
  // register the layer to the approach descriptor
  if (m_approachDescriptor) {
    approachDescriptor()->registerLayer(*this);
  }
}

std::shared_ptr<Acts::Layer>
Acts::DiscLayer::create(const variant_data& vardata)
{
  throw_assert(vardata.which() == 4, "Variant data must be map");
  const variant_map& data = boost::get<variant_map>(vardata);
  std::string        type = data.get<std::string>("type");
  throw_assert(type == "DiscLayer", "Type must be DiscLayer");

  variant_map payload = data.get<variant_map>("payload");

  auto trf = std::make_shared<const Transform3D>(
      from_variant<Transform3D>(payload.get<variant_map>("transform")));

  LayerType   laytyp;
  std::string laytyp_str = payload.get<std::string>("layer_type");
  if (laytyp_str == "active") {
    laytyp = active;
  } else if (laytyp_str == "passive") {
    laytyp = passive;
  } else { /*laytyp_str == "navigation"*/
    laytyp = navigation;
  }

  double thickness = payload.get<double>("thickness");
  double minR      = payload.get<double>("minR");
  double maxR      = payload.get<double>("maxR");
  double R         = (minR + maxR) / 2.;

  auto rbounds = std::make_shared<const RadialBounds>(minR, maxR);

  std::unique_ptr<SurfaceArray> sArray = nullptr;

  // only attempt to reover surface array if present
  if (payload.count("surfacearray") != 0u) {

    // get surface array transform
    const Transform3D& sa_trf = from_variant<Transform3D>(
        payload.get<variant_map>("surfacearray_transform"));
    const Transform3D& sa_itrf = sa_trf.inverse();

    // we need to reproduce the coordinate conversions
    auto g2l = [sa_trf](const Vector3D& pos) -> Vector2D {
      // @TODO: Check if - is right here, might be the other way round
      Vector3D loc = sa_trf * pos;
      return Vector2D(perp(loc), phi(loc));
    };
    auto l2g = [sa_itrf, R](const Vector2D& loc) -> Vector3D {
      return sa_itrf
          * Vector3D(R * std::cos(loc[0]), R * std::sin(loc[0]), loc[1]);
    };

    sArray = std::make_unique<SurfaceArray>(
        payload.at("surfacearray"),
        g2l,
        l2g,
        std::make_shared<const Transform3D>(sa_trf));
  }

  // @TODO: Implement ApproachDescriptor serialization
  return MutableLayerPtr(new DiscLayer(trf,
                                       rbounds,
                                       std::move(sArray),
                                       thickness,
                                       nullptr,  // std::move(ad),
                                       laytyp));
}

const Acts::DiscSurface&
Acts::DiscLayer::surfaceRepresentation() const
{
  return (*this);
}

Acts::DiscSurface&
Acts::DiscLayer::surfaceRepresentation()
{
  return (*this);
}

void
Acts::DiscLayer::buildApproachDescriptor()
{
  // delete it
  m_approachDescriptor.reset(nullptr);
  // take the boundary surfaces of the representving volume if they exist
  if (m_representingVolume != nullptr) {
    // get the boundary surfaces
    const std::vector<std::shared_ptr<const BoundarySurfaceT<AbstractVolume>>>&
        bSurfaces
        = m_representingVolume->boundarySurfaces();
    // fill in the surfaces into the vector
    std::vector<std::shared_ptr<const Surface>> aSurfaces;
    aSurfaces.push_back(
        bSurfaces.at(negativeFaceXY)->surfaceRepresentation().getSharedPtr());
    aSurfaces.push_back(
        bSurfaces.at(positiveFaceXY)->surfaceRepresentation().getSharedPtr());
    // create an ApproachDescriptor with Boundary surfaces
    m_approachDescriptor = std::make_unique<const GenericApproachDescriptor>(
        std::move(aSurfaces));
  } else {
    // create the new surfaces - positions first
    Vector3D aspPosition(center()
                         + 0.5 * thickness() * Surface::normal(center()));
    Vector3D asnPosition(center()
                         - 0.5 * thickness() * Surface::normal(center()));
    auto asnTransform
        = std::make_shared<const Transform3D>(Translation3D(asnPosition));
    auto aspTransform
        = std::make_shared<const Transform3D>(Translation3D(aspPosition));
    // create the vector
    std::vector<std::shared_ptr<const Surface>> aSurfaces;
    aSurfaces.push_back(
        Surface::makeShared<DiscSurface>(asnTransform, m_bounds));
    aSurfaces.push_back(
        Surface::makeShared<DiscSurface>(aspTransform, m_bounds));
    // create an ApproachDescriptor with standard surfaces surfaces - these
    // will be deleted by the approach descriptor
    m_approachDescriptor = std::make_unique<const GenericApproachDescriptor>(
        std::move(aSurfaces));
  }
  // @todo check if we can give the layer at curface creation
  for (auto& sfPtr : (m_approachDescriptor->containedSurfaces())) {
    if (sfPtr != nullptr) {
      auto& mutableSf = *(const_cast<Surface*>(sfPtr));
      mutableSf.associateLayer(*this);
    }
  }
}

Acts::variant_data
Acts::DiscLayer::toVariantData() const
{
  using namespace std::string_literals;
  variant_map payload;

  if (m_transform) {
    payload["transform"] = to_variant(*m_transform);
  }

  // we need to recover the bounds
  const AbstractVolume* absVol = representingVolume();
  throw_assert(absVol,
               "Cannot serialize DiscLayer without representing volume");
  auto cvBounds
      = dynamic_cast<const CylinderVolumeBounds*>(&absVol->volumeBounds());

  payload["minR"]      = cvBounds->innerRadius();
  payload["maxR"]      = cvBounds->outerRadius();
  payload["thickness"] = thickness();

  if (layerType() == active) {
    payload["layer_type"] = "active"s;
  } else if (layerType() == passive) {
    payload["layer_type"] = "passive"s;
  } else { /*layerType() == navigation*/
    payload["layer_type"] = "navigation"s;
  }

  if (m_surfaceArray) {
    payload["surfacearray"]           = m_surfaceArray->toVariantData();
    payload["surfacearray_transform"] = to_variant(m_surfaceArray->transform());
  }

  variant_map data;
  data["type"]    = "DiscLayer"s;
  data["payload"] = payload;
  return data;
}

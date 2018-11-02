// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderLayer.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Layers/CylinderLayer.hpp"
#include "Acts/Layers/GenericApproachDescriptor.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/InstanceFactory.hpp"
#include "Acts/Utilities/VariantData.hpp"
#include "Acts/Volumes/AbstractVolume.hpp"
#include "Acts/Volumes/BoundarySurfaceFace.hpp"
#include "Acts/Volumes/CylinderVolumeBounds.hpp"

using Acts::VectorHelpers::phi;

Acts::CylinderLayer::CylinderLayer(
    const std::shared_ptr<const Transform3D>&    transform,
    const std::shared_ptr<const CylinderBounds>& cBounds,
    std::unique_ptr<SurfaceArray>                surfaceArray,
    double                                       thickness,
    std::unique_ptr<ApproachDescriptor>          ades,
    LayerType                                    laytyp)
  : CylinderSurface(transform, cBounds)
  , Layer(std::move(surfaceArray), thickness, std::move(ades), laytyp)
{
  // create the representing volume
  CylinderVolumeBounds* cvBounds
      = new CylinderVolumeBounds(cBounds->r() - 0.5 * thickness,
                                 cBounds->r() + 0.5 * thickness,
                                 cBounds->halflengthZ());
  m_representingVolume
      = new AbstractVolume(transform, VolumeBoundsPtr(cvBounds));
  // associate the layer to the surface
  CylinderSurface::associateLayer(*this);
  // an approach descriptor is automatically created if there's a surface array
  if (!m_approachDescriptor && m_surfaceArray) {
    buildApproachDescriptor();
  }
  // register the layer to the approach descriptor surfaces
  if (m_approachDescriptor) {
    approachDescriptor()->registerLayer(*this);
  }
}

std::shared_ptr<Acts::Layer>
Acts::CylinderLayer::create(const variant_data& vardata)
{
  throw_assert(vardata.which() == 4, "Variant data must be map");
  const variant_map& data = boost::get<variant_map>(vardata);
  std::string        type = data.get<std::string>("type");
  throw_assert(type == "CylinderLayer", "Type must be CylinderLayer");

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

  InstanceFactory    factory;
  const variant_map& var_bounds = payload.get<variant_map>("cylinder_bounds");

  auto cbounds = std::dynamic_pointer_cast<const CylinderBounds>(
      factory.surfaceBounds(var_bounds.get<std::string>("type"), var_bounds));

  std::unique_ptr<SurfaceArray> sArray = nullptr;

  // only attempt to reover surface array if present
  if (payload.count("surfacearray") != 0u) {

    // get surface array transform
    const Transform3D& sa_trf = from_variant<Transform3D>(
        payload.get<variant_map>("surfacearray_transform"));
    const Transform3D& sa_itrf = sa_trf.inverse();

    // we need to reproduce the coordinate conversions
    double R      = cbounds->r();
    double avgPhi = cbounds->averagePhi();
    auto g2l      = [sa_trf, avgPhi](const Vector3D& pos) -> Vector2D {
      // @TODO: Check if - is right here, might be the other way round
      Vector3D loc = sa_trf * pos;
      return Vector2D(phi(loc) - avgPhi, loc.z());
    };
    auto l2g = [sa_itrf, R, avgPhi](const Vector2D& loc) -> Vector3D {
      return sa_itrf * Vector3D(R * std::cos(loc[0] + avgPhi),
                                R * std::sin(loc[0] + avgPhi),
                                loc[1]);
    };

    sArray = std::make_unique<SurfaceArray>(
        payload.at("surfacearray"),
        g2l,
        l2g,
        std::make_shared<const Transform3D>(sa_trf));
  }

  // @TODO: Implement ApproachDescriptor serialization
  return MutableLayerPtr(new CylinderLayer(trf,
                                           cbounds,
                                           std::move(sArray),
                                           thickness,
                                           nullptr,  // std::move(ad),
                                           laytyp));
}

const Acts::CylinderSurface&
Acts::CylinderLayer::surfaceRepresentation() const
{
  return (*this);
}

Acts::CylinderSurface&
Acts::CylinderLayer::surfaceRepresentation()
{
  return (*this);
}

void
Acts::CylinderLayer::buildApproachDescriptor()
{
  // delete it
  m_approachDescriptor.reset(nullptr);
  // delete the surfaces
  // take the boundary surfaces of the representving volume if they exist
  if (m_representingVolume != nullptr) {
    // get the boundary surfaces
    const std::vector<std::shared_ptr<const BoundarySurfaceT<AbstractVolume>>>&
        bSurfaces
        = m_representingVolume->boundarySurfaces();

    // fill in the surfaces into the vector
    std::vector<std::shared_ptr<const Surface>> aSurfaces;
    if (bSurfaces.size() > size_t(tubeOuterCover)) {
      aSurfaces.push_back(
          bSurfaces.at(tubeInnerCover)->surfaceRepresentation().getSharedPtr());
    }
    aSurfaces.push_back(
        bSurfaces.at(tubeOuterCover)->surfaceRepresentation().getSharedPtr());
    // create an ApproachDescriptor with Boundary surfaces
    m_approachDescriptor = std::make_unique<const GenericApproachDescriptor>(
        std::move(aSurfaces));
  } else {
    // create the new surfaces
    std::vector<std::shared_ptr<const Acts::Surface>> aSurfaces;
    aSurfaces.push_back(
        Surface::makeShared<CylinderSurface>(m_transform,
                                             m_bounds->r() - 0.5 * thickness(),
                                             m_bounds->halflengthZ()));
    aSurfaces.push_back(
        Surface::makeShared<CylinderSurface>(m_transform,
                                             m_bounds->r() + 0.5 * thickness(),
                                             m_bounds->halflengthZ()));
    // create an ApproachDescriptor with standard surfaces surfaces - these will
    // be deleted by the approach descriptor
    m_approachDescriptor = std::make_unique<const GenericApproachDescriptor>(
        std::move(aSurfaces));
  }
  for (auto& sfPtr : (m_approachDescriptor->containedSurfaces())) {
    if (sfPtr != nullptr) {
      auto& mutableSf = *(const_cast<Surface*>(sfPtr));
      mutableSf.associateLayer(*this);
    }
  }
}

Acts::variant_data
Acts::CylinderLayer::toVariantData() const
{
  using namespace std::string_literals;
  variant_map payload;

  if (m_transform == nullptr) {
    payload["transform"] = to_variant(Transform3D::Identity());
  } else {
    payload["transform"] = to_variant(*m_transform);
  }

  // we need to recover the bounds
  const AbstractVolume* absVol = representingVolume();
  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  auto cvBounds
      = dynamic_cast<const CylinderVolumeBounds*>(&absVol->volumeBounds());

  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  double cylR = cvBounds->innerRadius() + 0.5 * thickness();
  double hlZ  = cvBounds->halflengthZ();

  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  CylinderBounds cylBounds(cylR, hlZ);

  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  const variant_data bounds = cylBounds.toVariantData();
  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  payload["cylinder_bounds"] = bounds;

  payload["thickness"] = thickness();

  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
  if (layerType() == active) {
    payload["layer_type"] = "active"s;
  } else if (layerType() == passive) {
    payload["layer_type"] = "passive"s;
  } else { /*layerType() == navigation*/
    payload["layer_type"] = "navigation"s;
  }

  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

  // this lacks localToGlobal and globalToLocal
  if (m_surfaceArray) {
    payload["surfacearray"]           = m_surfaceArray->toVariantData();
    payload["surfacearray_transform"] = to_variant(m_surfaceArray->transform());
  }

  variant_map data;
  data["type"]    = "CylinderLayer"s;
  data["payload"] = payload;
  return data;
}

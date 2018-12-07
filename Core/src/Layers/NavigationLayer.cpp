// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// NavigationLayer.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Layers/NavigationLayer.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/InstanceFactory.hpp"
#include "Acts/Utilities/VariantData.hpp"

Acts::NavigationLayer::NavigationLayer(
    std::shared_ptr<const Surface> surfaceRepresentation,
    double                         thickness)
  : Acts::Layer(nullptr)
  , m_surfaceRepresentation(std::move(surfaceRepresentation))
{
  m_layerThickness = thickness;
  m_layerType      = navigation;
}

std::shared_ptr<const Acts::Layer>
Acts::NavigationLayer::create(const variant_data& vardata)
{
  throw_assert(vardata.which() == 4, "Variant data must be map");
  const variant_map& data = boost::get<variant_map>(vardata);
  std::string        type = data.get<std::string>("type");
  throw_assert(type == "NavigationLayer", "Type must be NavigationLayer");

  const variant_map& payload = data.get<variant_map>("payload");

  double          thickness = payload.get<double>("thickness");
  InstanceFactory factory;

  const variant_map& var_surface
      = payload.get<variant_map>("surface_representation");
  std::shared_ptr<const Surface> surface_ptr
      = factory.surface(var_surface.get<std::string>("type"), var_surface);

  return NavigationLayer::create(std::move(surface_ptr), thickness);
}

Acts::NavigationLayer::~NavigationLayer() = default;

Acts::variant_data
Acts::NavigationLayer::toVariantData() const
{
  using namespace std::string_literals;
  variant_map payload;

  payload["surface_representation"] = m_surfaceRepresentation->toVariantData();
  payload["thickness"]              = thickness();

  variant_map data;
  data["type"]    = "NavigationLayer"s;
  data["payload"] = payload;
  return data;
}

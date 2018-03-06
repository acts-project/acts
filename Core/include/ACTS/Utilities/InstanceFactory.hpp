// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_PLUGINS_JSONPLUGIN_INSTANCEFACTORY_H
#define ACTS_PLUGINS_JSONPLUGIN_INSTANCEFACTORY_H 1

#include <boost/functional/value_factory.hpp>
#include <map>
#include "ACTS/Utilities/VariantData.hpp"
#include "ACTS/Utilities/ThrowAssert.hpp"

#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Surfaces/TriangleBounds.hpp"
#include "ACTS/Surfaces/DiamondBounds.hpp"
#include "ACTS/Surfaces/EllipseBounds.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/TrapezoidBounds.hpp"

namespace Acts {

using SurfaceBoundsPtr     = std::shared_ptr<const SurfaceBounds>;
using PlanarBoundsPtr     = std::shared_ptr<const PlanarBounds>;
using SurfaceBoundsFactory = std::function<SurfaceBoundsPtr(const variant_data&)>;

class InstanceFactory
{
public:
  InstanceFactory()
  {
    // set up map to store factories
    m_surfaceBounds["RectangleBounds"] = [](auto const& data) {
      return std::make_shared<const RectangleBounds>(data);
    };
    m_surfaceBounds["TriangleBounds"] = [](auto const& data) {
      return std::make_shared<const TriangleBounds>(data);
    };
    m_surfaceBounds["DiamondBounds"] = [](auto const& data) {
      return std::make_shared<const DiamondBounds>(data);
    };
    m_surfaceBounds["EllipseBounds"] = [](auto const& data) {
      return std::make_shared<const EllipseBounds>(data);
    };
    m_surfaceBounds["TrapezoidBounds"] = [](auto const& data) {
      return std::make_shared<const TrapezoidBounds>(data);
    };
    m_surfaceBounds["CylinderBounds"] = [](auto const& data) {
      return std::make_shared<const CylinderBounds>(data);
    };

  }

  PlanarBoundsPtr
  planarBounds(const std::string& cname, const variant_data& data)
  {
    SurfaceBoundsPtr srfBnd_ptr = surfaceBounds(cname, data);
    PlanarBoundsPtr plnBnd_ptr = std::dynamic_pointer_cast<const PlanarBounds>(srfBnd_ptr);
    throw_assert(plnBnd_ptr, "Conversion to PlanarBounds failed");
    return plnBnd_ptr;
  }

  SurfaceBoundsPtr
  surfaceBounds(const std::string& cname, const variant_data& data)
  {
    throw_assert(m_surfaceBounds.count(cname), "No factory found for class "+cname);
    return m_surfaceBounds.at(cname)(data);
  }

private:
  std::map<std::string, SurfaceBoundsFactory> m_surfaceBounds;
};

}  // namespace Acts

#endif

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

namespace Acts {

using PlanarBoundsPtr     = std::shared_ptr<const PlanarBounds>;
using PlanarBoundsFactory = std::function<PlanarBoundsPtr(const variant_data&)>;

class InstanceFactory
{
public:
  InstanceFactory()
  {
    // set up map to store factories
    m_planarBounds["RectangleBounds"] = [](auto const& data) {
      return std::make_shared<const RectangleBounds>(data);
    };
    m_planarBounds["TriangleBounds"] = [](auto const& data) {
      return std::make_shared<const TriangleBounds>(data);
    };
  }

  PlanarBoundsPtr
  planarBounds(const std::string& cname, const variant_data& data)
  {
    throw_assert(m_planarBounds.count(cname), "No factory found for class "+cname);
    return m_planarBounds.at(cname)(data);
  }

private:
  std::map<std::string, PlanarBoundsFactory> m_planarBounds;
};

}  // namespace Acts

#endif

// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <boost/functional/value_factory.hpp>
#include <map>
#include <memory>
#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/VariantData.hpp"

#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/DiscTrapezoidalBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Surfaces/TriangleBounds.hpp"

namespace Acts {

/// Class that creates instances for certain types of objects
/// from @c variant_data input. This is used to avoid duplicating
/// this type of code anywhere we construct from @c variant_data.
class InstanceFactory
{
  // internal typedefs
  using SurfaceBoundsPtr = std::shared_ptr<const SurfaceBounds>;
  using PlanarBoundsPtr  = std::shared_ptr<const PlanarBounds>;
  using SurfaceBoundsFactory
      = std::function<SurfaceBoundsPtr(const variant_data&)>;
  using SurfaceFactory
      = std::function<std::shared_ptr<Surface>(const variant_data&)>;

public:
  /// Default constructor. Sets up a map to lambdas which return
  /// pointers to the newly constructor instances.
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
    m_surfaceBounds["RadialBounds"] = [](auto const& data) {
      return std::make_shared<const RadialBounds>(data);
    };
    m_surfaceBounds["DiscTrapezoidalBounds"] = [](auto const& data) {
      return std::make_shared<const DiscTrapezoidalBounds>(data);
    };

    m_surfaces["PlaneSurface"] = [](auto const& data) {
      return Surface::makeShared<PlaneSurface>(data);
    };
  }

  /// Method to produce surface bounds type objects
  /// @param cname The class name
  /// @param data The @c variant_data to construct from
  /// @return Shared pointer to the created surfacebounds
  SurfaceBoundsPtr
  surfaceBounds(const std::string& cname, const variant_data& data) const
  {
    throw_assert(m_surfaceBounds.count(cname),
                 "No factory found for class " + cname);
    return m_surfaceBounds.at(cname)(data);
  }

  /// Method to produce planar bounds type objects
  /// @param cname The class name
  /// @param data The @c variant_data to construct from
  /// @return Shared pointer to the created planar bounds
  PlanarBoundsPtr
  planarBounds(const std::string& cname, const variant_data& data) const
  {
    SurfaceBoundsPtr srfBnd_ptr = surfaceBounds(cname, data);
    PlanarBoundsPtr  plnBnd_ptr
        = std::dynamic_pointer_cast<const PlanarBounds>(srfBnd_ptr);
    throw_assert(plnBnd_ptr,
                 "Conversion to PlanarBounds failed for name " + cname);
    return plnBnd_ptr;
  }

  /// Method to produce disc bounds type objects
  /// @param cname The class name
  /// @param data The @c variant_data to construct from
  /// @return Shared pointer to the created disc bounds
  std::shared_ptr<const DiscBounds>
  discBounds(const std::string& cname, const variant_data& data) const
  {
    SurfaceBoundsPtr                  srfBnd_ptr = surfaceBounds(cname, data);
    std::shared_ptr<const DiscBounds> discBnd_ptr
        = std::dynamic_pointer_cast<const DiscBounds>(srfBnd_ptr);
    throw_assert(discBnd_ptr,
                 "Conversion to DiscBounds failed for name " + cname);
    return discBnd_ptr;
  }

  /// Method to produce surface objects
  /// @param cname The class name
  /// @param data The @c variant_data to construct from
  /// @return Pointer to the created surface
  std::shared_ptr<Surface>
  surface(const std::string& cname, const variant_data& data) const
  {
    throw_assert(m_surfaces.count(cname),
                 "No factory found for class " + cname);
    return m_surfaces.at(cname)(data);
  }

private:
  std::map<std::string, SurfaceBoundsFactory> m_surfaceBounds;
  std::map<std::string, SurfaceFactory>       m_surfaces;
};

}  // namespace Acts
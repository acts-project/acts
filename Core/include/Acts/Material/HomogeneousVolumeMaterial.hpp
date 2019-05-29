// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// HomogeneousVolumeMaterial.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class HomogeneousVolumeMaterial
///
/// It extends the IVolumeMaterial base class to describe a simple
/// homogeneous material in a volume
class HomogeneousVolumeMaterial : public IVolumeMaterial {
 public:
  /// Default Constructor - defaulted
  HomogeneousVolumeMaterial() = default;

  /// Explicit constructor
  ///
  /// @param material is the material held by this
  HomogeneousVolumeMaterial(const Material& material) : m_material(material) {}

  /// Access to actual material
  ///
  /// @param position is the request position for the material call
  /// @todo interface to change including 'cell'
  const Material& material(const Vector3D& /*position*/) const final {
    return m_material;
  }

 private:
  Material m_material;
};

}  // namespace Acts

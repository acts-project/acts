// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
  HomogeneousVolumeMaterial(const Material& material);

  /// Copy Constructor
  ///
  /// @param hvm is the source material
  HomogeneousVolumeMaterial(const HomogeneousVolumeMaterial& hvm) = default;

  /// Copy Move Constructor
  ///
  /// @param hvm is the source material
  HomogeneousVolumeMaterial(HomogeneousVolumeMaterial&& hvm) = default;

  /// Destructor
  ~HomogeneousVolumeMaterial() override = default;

  /// Assignment operator
  ///
  /// @param hvm is the source material
  HomogeneousVolumeMaterial& operator=(const HomogeneousVolumeMaterial& hvm) =
      default;

  /// Equality operator
  ///
  /// @param hvm is the source material
  bool operator==(const HomogeneousVolumeMaterial& hvm) const;

  /// Access to actual material
  ///
  /// @param position is the request position for the material call
  /// @todo interface to change including 'cell'
  const Material material(const Vector3D& /*position*/) const final;

  /// Output Method for std::ostream
  ///
  /// @param sl The outoput stream
  std::ostream& toStream(std::ostream& sl) const;

 private:
  Material m_material = Material();
};

inline const Material HomogeneousVolumeMaterial::material(
    const Vector3D& /*position*/) const {
  return (m_material);
}

inline bool HomogeneousVolumeMaterial::operator==(
    const HomogeneousVolumeMaterial& hvm) const {
  return (m_material == hvm.m_material);
}

}  // namespace Acts

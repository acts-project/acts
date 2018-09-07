// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GridSurfaceMaterial.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <vector>
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class GridSurfaceMaterial
///
/// It extends the SurfaceMaterial base class and describes a simple homogeneous
/// material descriptions

class GridSurfaceMaterial : public SurfaceMaterial
{
public:
  /// Default Constructor - defaulted
  GridSurfaceMaterial() = default;

  /// Explicit constructor
  ///
  /// @param full are the full material properties
  /// @param splitFactor is the split for pre/post update
  GridSurfaceMaterial(const MaterialProperties& full, double splitFactor = 1.);

  /// Copy Constructor
  ///
  /// @param hsm is the source material
  GridSurfaceMaterial(const GridSurfaceMaterial& hsm) = default;

  /// Copy Move Constructor
  ///
  /// @param hsm is the source material
  GridSurfaceMaterial(GridSurfaceMaterial&& hsm) = default;

  /// Destructor
  ~GridSurfaceMaterial() override = default;

  /// Assignment operator
  ///
  /// @param hsm is the source material
  GridSurfaceMaterial&
  operator=(const GridSurfaceMaterial& hsm)
      = default;

  /// Assignment Move operator
  ///
  /// @param hsm is the source material
  GridSurfaceMaterial&
  operator=(GridSurfaceMaterial&& hsm)
      = default;

  /// Scale operator
  /// - it is effectively a thickness scaling
  ///
  /// @param scale is the scale factor
  GridSurfaceMaterial&
  operator*=(double scale) final;

  /// Equality operator
  ///
  /// @param hsm is the source material
  bool
  operator==(const GridSurfaceMaterial& hsm) const;

  /// @copydoc SurfaceMaterial::materialProperties(const Vector2D&)
  ///
  /// @note the input parameter is ignored
  const MaterialProperties&
  materialProperties(const Vector2D& lp) const final;

  /// @copydoc SurfaceMaterial::materialProperties(const Vector3D&)
  ///
  /// @note the input parameter is ignored
  const MaterialProperties&
  materialProperties(const Vector3D& gp) const final;

  /// @copydoc SurfaceMaterial::materialProperties(size_t, size_t)
  ///
  /// @param ib0 The bin at local 0 for retrieving the material
  /// @param ib1 The bin at local 1 for retrieving the material
  ///
  /// @note the input parameter is ignored
  const MaterialProperties&
  materialProperties(size_t ib0, size_t ib1) const final;

  /// The inherited methods - for materialProperties access
  using SurfaceMaterial::materialProperties;

  /// The interited methods - for scale access
  using SurfaceMaterial::factor;

  /// Output Method for std::ostream
  ///
  /// @param sl The outoput stream
  std::ostream&
  dump(std::ostream& sl) const final;

private:
  /// The five different MaterialProperties
  MaterialProperties m_fullMaterial = MaterialProperties();
};

inline const MaterialProperties&
GridSurfaceMaterial::materialProperties(const Vector2D& /*lp*/) const
{
  return (m_fullMaterial);
}

inline const MaterialProperties&
GridSurfaceMaterial::materialProperties(const Vector3D& /*gp*/) const
{
  return (m_fullMaterial);
}

inline const MaterialProperties&
    GridSurfaceMaterial::materialProperties(size_t /*ib0*/,
                                            size_t /*ib1*/) const
{
  return (m_fullMaterial);
}

inline bool
GridSurfaceMaterial::operator==(const GridSurfaceMaterial& hsm) const
{
  return (m_fullMaterial == hsm.m_fullMaterial);
}

}  // namespace
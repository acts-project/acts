// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"

#include <iosfwd>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Acts {

/// @addtogroup material
/// @{

///
/// @brief Surface material placeholder carrying a two-dimensional axis spec
///
/// The ProtoSurfaceMaterial class acts as a proxy to the SurfaceMaterial
/// to mark the layers and surfaces on which the material should be mapped on
/// at construction time of the geometry and to hand over the granularity of
/// the material map with a MultiAxisSpec2D. Deferred axes are resolved
/// against the surface during material mapping. A single bin in both
/// directions requests homogeneous material.
class ProtoSurfaceMaterial final : public ISurfaceMaterial {
 public:
  /// Construct a homogeneous placeholder with one deferred bin per axis
  ProtoSurfaceMaterial() = default;

  /// Constructor with MultiAxisSpec2D
  /// @param binning a binning description for the material map binning
  /// @param materialKey Optional stable identity of the material assignment
  /// @param mappingType is the type of surface mapping associated to the surface
  explicit ProtoSurfaceMaterial(
      const MultiAxisSpec2D& binning,
      MappingType mappingType = MappingType::Default,
      std::optional<std::string> materialKey = std::nullopt)
      : ISurfaceMaterial(1., mappingType),
        m_binning(binning),
        m_materialKey(std::move(materialKey)) {
    if (m_materialKey && m_materialKey->empty()) {
      throw std::invalid_argument("Material key must not be empty");
    }
  }

  /// Copy constructor
  ///
  /// @param smproxy The source proxy
  ProtoSurfaceMaterial(const ProtoSurfaceMaterial& smproxy) = default;

  /// Copy move constructor
  ///
  /// @param smproxy The source proxy
  ProtoSurfaceMaterial(ProtoSurfaceMaterial&& smproxy) noexcept = default;

  /// Destructor
  /// Defined out of line so ActsCore owns the vtable and RTTI used when
  /// dispatching material across shared-library boundaries.
  ~ProtoSurfaceMaterial() override;

  /// Assignment operator
  ///
  /// @param smproxy The source proxy
  /// @return Reference to this object
  ProtoSurfaceMaterial& operator=(const ProtoSurfaceMaterial& smproxy) =
      default;

  /// Assignment move operator
  ///
  /// @param smproxy The source proxy
  /// @return Reference to this object
  ProtoSurfaceMaterial& operator=(ProtoSurfaceMaterial&& smproxy) noexcept =
      default;

  /// Scale operation - dummy implementation
  ///
  /// @return Reference to this object
  ProtoSurfaceMaterial& scale(double /*factor*/) override { return *this; }

  /// Stable identity of the material assignment, if configured
  /// @return Optional stable material assignment key
  const std::optional<std::string>& materialKey() const {
    return m_materialKey;
  }

  /// Return the two-dimensional binning specification
  /// @return Reference to the binning
  const MultiAxisSpec2D& binning() const { return m_binning; }

  /// Return method for full material description of the Surface - from local
  /// coordinates
  ///
  /// @return will return dummy material
  const MaterialSlab& materialSlab(const Vector2& /*lp*/) const override {
    return m_materialSlab;
  }

  /// @copydoc ISurfaceMaterial::localAxisDirections() const
  /// A placeholder has no lookup grid. Axis directions are validated and
  /// ordered when the binning is resolved for material mapping.
  std::vector<AxisDirection> localAxisDirections() const override { return {}; }

  using ISurfaceMaterial::materialSlab;

  /// Output Method for std::ostream, to be overloaded by child classes
  ///
  /// @param sl is the output stream
  /// @return The output stream
  std::ostream& toStream(std::ostream& sl) const override {
    sl << "Acts::ProtoSurfaceMaterial : " << std::endl;
    sl << m_binning << std::endl;
    return sl;
  }

 private:
  /// A binning description
  MultiAxisSpec2D m_binning{
      {AxisSpec::DeferredEquidistant(1), AxisSpec::DeferredEquidistant(1)}};

  std::optional<std::string> m_materialKey;

  /// Dummy material properties
  MaterialSlab m_materialSlab = MaterialSlab::Nothing();
};

/// @}

}  // namespace Acts

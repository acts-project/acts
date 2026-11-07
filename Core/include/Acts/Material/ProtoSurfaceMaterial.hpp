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

#include <array>
#include <iosfwd>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace Acts {

/// @addtogroup material
/// @{

///
/// @brief Proxy to surface material carrying a two-dimensional binning spec
///
/// The ProtoSurfaceMaterial class acts as a proxy to the SurfaceMaterial
/// to mark the layers and surfaces on which the material should be mapped on
/// at construction time of the geometry and to hand over the granularity
/// of the material map. Deferred axes are resolved against the surface during
/// material mapping; one bin in each direction describes homogeneous material.
class ProtoSurfaceMaterial : public ISurfaceMaterial {
 public:
  /// Constructor without a binning spec - homogeneous material
  ProtoSurfaceMaterial() = default;

  /// Constructor with a two-dimensional binning spec
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
  ~ProtoSurfaceMaterial() override = default;

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
  ProtoSurfaceMaterial& scale(double /*factor*/) final { return (*this); }

  /// Stable identity of the material assignment, if configured
  /// @return Optional stable material assignment key
  const std::optional<std::string>& materialKey() const {
    return m_materialKey;
  }

  /// Return the multi-axis binning spec
  /// @return Reference to the binning
  const MultiAxisSpec2D& binning() const { return m_binning; }

  /// Return method for full material description of the Surface - from local
  /// coordinates
  ///
  /// @return will return dummy material
  const MaterialSlab& materialSlab(const Vector2& /*lp*/) const final {
    return (m_materialSlab);
  }

  /// @copydoc ISurfaceMaterial::localAxisDirections() const
  std::vector<AxisDirection> localAxisDirections() const final { return {}; }

  using ISurfaceMaterial::materialSlab;

  /// Output Method for std::ostream, to be overloaded by child classes
  ///
  /// @param sl is the output stream
  /// @return The output stream
  std::ostream& toStream(std::ostream& sl) const final {
    sl << "Acts::ProtoSurfaceMaterial : " << std::endl;
    sl << m_binning << std::endl;
    return sl;
  }

 private:
  /// A binning description
  MultiAxisSpec2D m_binning{std::array{AxisSpec::DeferredEquidistant(1),
                                       AxisSpec::DeferredEquidistant(1)}};

  std::optional<std::string> m_materialKey;

  /// Dummy material properties
  MaterialSlab m_materialSlab = MaterialSlab::Nothing();
};

/// @}

}  // namespace Acts

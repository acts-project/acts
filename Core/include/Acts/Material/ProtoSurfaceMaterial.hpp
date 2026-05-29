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
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

#include <array>
#include <iosfwd>
#include <vector>

namespace Acts {

/// @addtogroup material
/// @{

/// Convert a BinUtility to a pair of two DirectedProtoAxis.
/// For a 1D BinUtility a single-bin dummy axis is appended as the second axis.
/// This helper is intended for migration code that still works with BinUtility.
std::array<DirectedProtoAxis, 2> protoAxesFromBinUtility(const BinUtility& bu);

///
/// @brief Prototype surface material described by two @ref DirectedProtoAxis
///        objects.
///
/// The ProtoSurfaceMaterial class acts as a proxy to the SurfaceMaterial
/// to mark surfaces on which the material should be mapped at geometry
/// construction time. It carries the binning description with two directed
/// proto axes. A conceptually 1D surface uses a second axis with a single bin.
///
class ProtoSurfaceMaterial : public ISurfaceMaterial {
 public:
  /// Default constructor — creates a homogeneous (single-bin) placeholder.
  ProtoSurfaceMaterial();

  /// Constructor with two directed proto axes.
  ///
  /// @param axes the 2D binning description
  /// @param mappingType is the type of surface mapping associated to the surface
  explicit ProtoSurfaceMaterial(std::array<DirectedProtoAxis, 2> axes,
                                MappingType mappingType = MappingType::Default);

  /// @deprecated Use the std::array<DirectedProtoAxis, 2> constructor instead.
  ///
  /// Constructor from a BinUtility.
  ///
  /// @param binUtil the binning utility (converted to two DirectedProtoAxis)
  /// @param mappingType is the type of surface mapping associated to the surface
  [[deprecated(
      "Use ProtoSurfaceMaterial(std::array<DirectedProtoAxis, 2>) "
      "instead")]]
  explicit ProtoSurfaceMaterial(const BinUtility& binUtil,
                                MappingType mappingType = MappingType::Default);

  /// @deprecated Use the std::array<DirectedProtoAxis, 2> constructor instead.
  ///
  /// Constructor from a vector of DirectedProtoAxis objects (was the
  /// ProtoGridSurfaceMaterial path). Only the first two axes are used; if only
  /// one is provided a single-bin dummy axis is added.
  ///
  /// @param axes vector of directed proto axes
  /// @param mappingType is the type of surface mapping associated to the surface
  [[deprecated(
      "Use ProtoSurfaceMaterial(std::array<DirectedProtoAxis, 2>) "
      "instead")]]
  explicit ProtoSurfaceMaterial(const std::vector<DirectedProtoAxis>& axes,
                                MappingType mappingType = MappingType::Default);

  /// Copy / move / destructor — all defaulted.
  ProtoSurfaceMaterial(const ProtoSurfaceMaterial&) = default;
  ProtoSurfaceMaterial(ProtoSurfaceMaterial&&) noexcept = default;
  ~ProtoSurfaceMaterial() override = default;
  ProtoSurfaceMaterial& operator=(const ProtoSurfaceMaterial&) = default;
  ProtoSurfaceMaterial& operator=(ProtoSurfaceMaterial&&) noexcept = default;

  /// Scale operation — no-op for proto material.
  ProtoSurfaceMaterial& scale(double /*factor*/) final { return *this; }

  /// Return the two directed proto axes.
  /// @return const reference to the axes array
  const std::array<DirectedProtoAxis, 2>& axes() const { return m_axes; }

  /// @deprecated Use axes() instead.
  ///
  /// Reconstruct a BinUtility from the internal axes.
  /// Single-bin dummy axes are omitted so the result matches what was
  /// originally provided to the deprecated BinUtility constructor.
  ///
  /// @return reconstructed BinUtility
  [[deprecated("Use axes() instead")]] BinUtility binning() const;

  /// Return method for full material description — returns dummy material.
  const MaterialSlab& materialSlab(const Vector2& /*lp*/) const final {
    return m_materialSlab;
  }

  /// @copydoc ISurfaceMaterial::localAxisDirections() const
  std::vector<AxisDirection> localAxisDirections() const final;

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector3&) const
  [[deprecated(
      "Use materialSlab(const Vector2& lp) with a prior "
      "Surface::globalToLocal() call instead")]] const MaterialSlab&
  materialSlab(const Vector3& /*gp*/) const final {
    return m_materialSlab;
  }

  using ISurfaceMaterial::materialSlab;

  /// Output Method for std::ostream.
  std::ostream& toStream(std::ostream& sl) const final;

 private:
  std::array<DirectedProtoAxis, 2> m_axes;
  MaterialSlab m_materialSlab = MaterialSlab::Nothing();
};

/// @deprecated Use ProtoSurfaceMaterial instead. Will be removed in a future release.
using ProtoGridSurfaceMaterial
    [[deprecated("Use ProtoSurfaceMaterial instead")]] = ProtoSurfaceMaterial;

/// @}

}  // namespace Acts

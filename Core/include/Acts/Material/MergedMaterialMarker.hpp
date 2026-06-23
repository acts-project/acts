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

#include <iosfwd>
#include <vector>

namespace Acts {

/// @ingroup material
///
/// Sentinel surface material used by `Portal::merge` in its "keep going" mode.
///
/// When two portal surfaces that both carry material (or one of them does) are
/// merged, the original material cannot be transferred onto the (larger) merged
/// surface. Instead of aborting the construction, the merge can be configured
/// to discard the input material and tag the merged surface with this marker.
///
/// The marker carries no physical material -- it always returns
/// @ref MaterialSlab::Nothing() -- but its presence makes the lossy merge
/// discoverable downstream (e.g. when inspecting or writing out the geometry).
class MergedMaterialMarker final : public ISurfaceMaterial {
 public:
  /// Default constructor
  MergedMaterialMarker() = default;

  /// Destructor
  ~MergedMaterialMarker() override = default;

  /// Scale operator -- no-op, the marker carries no material
  /// @param factor is the scale factor (ignored)
  /// @return Reference to this marker
  MergedMaterialMarker& scale(double factor) override;

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector2&) const
  ///
  /// @note the input parameter is ignored, always returns
  ///       @ref MaterialSlab::Nothing()
  const MaterialSlab& materialSlab(const Vector2& lp = Vector2{
                                       0., 0.}) const override;

  /// @copydoc ISurfaceMaterial::localAxisDirections() const
  std::vector<AxisDirection> localAxisDirections() const override;

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector3&) const
  ///
  /// @note the input parameter is ignored, always returns
  ///       @ref MaterialSlab::Nothing()
  [[deprecated(
      "Use materialSlab(const Vector2& lp) with a prior "
      "Surface::globalToLocal() call instead")]] const MaterialSlab&
  materialSlab(const Vector3& gp) const override;

  // Inherit additional materialSlab overloads from base class
  using ISurfaceMaterial::materialSlab;

  // Inherit the scale-access helper from the base class
  using ISurfaceMaterial::factor;

  /// Output Method for std::ostream
  /// @param sl The output stream
  /// @return Reference to the output stream for chaining
  std::ostream& toStream(std::ostream& sl) const override;

 private:
  /// The marker carries no material
  MaterialSlab m_slab = MaterialSlab::Nothing();
};

}  // namespace Acts

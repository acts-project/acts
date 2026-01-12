// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/Material.hpp"

#include <iosfwd>
#include <limits>
#include <utility>
#include <vector>

// Tell the compiler to optimize the containing block assuming that
// FP may trap.  This is sometimes needed with clang to avoid spurious FPEs
// resulting from auto-vectorization.

#if defined(__clang__) && defined(__x86_64__)
#pragma float_control(push)
#pragma float_control(except, on)
#endif

namespace Acts {

/// Material description for an object with defined thickness.
///
/// @ingroup material
///
/// This is intended to describe concrete surface materials.
///
/// @see Material for a description of the available parameters.
class MaterialSlab {
 public:
  /// Create a material slab with no material content
  /// @return Empty material slab with zero thickness and no material
  static constexpr MaterialSlab Nothing() {
    return MaterialSlab(Material::Vacuum(), 0, false);
  }

  /// Create a vacuum material slab with specified thickness
  /// @param thickness The thickness of the vacuum region
  /// @return Vacuum material slab with the given thickness
  static constexpr MaterialSlab Vacuum(float thickness) {
    return MaterialSlab(Material::Vacuum(), thickness, false);
  }

  /// Combine material properties of two layers by averaging them.
  ///
  /// @param layerA Input layer A to average over.
  /// @param layerB Input layer B to average over.
  ///
  /// @return The resulting object has the combined thickness of all layers but just
  ///         one set of appropriately averaged material constants.
  static MaterialSlab combineLayers(const MaterialSlab& layerA,
                                    const MaterialSlab& layerB);

  /// Compute the average properties for a combined slab of two materials.
  ///
  /// The averaged material slab has the combined thickness of the two input
  /// slabs and assumes the two input materials are homogeneously and
  /// continuously mixed throughout the slab.
  ///
  /// @param slab1 Properties of the first material slab
  /// @param material2 Properties of the second material
  /// @param thickness2 Thickness of the second material slab. Can be negative to
  ///                   subtract the second material from the first slab.
  ///
  /// @returns Material slab with the combined thickness and average parameters
  static MaterialSlab combine(const MaterialSlab& slab1,
                              const Material& material2, float thickness2);

  /// Combine material properties of multiple layers by averaging them.
  ///
  /// @param layers Input layers to average over.
  ///
  /// @return The resulting object has the combined thickness of all layers but just
  ///         one set of appropriately averaged material constants.
  static MaterialSlab combineLayers(const std::vector<MaterialSlab>& layers);

  /// Default constructor.
  ///
  /// TODO consider removing. currently needed for default construction in grids
  constexpr MaterialSlab() : m_material(Material::Vacuum()) {}

  /// Construct from material description.
  ///
  /// @param material  is the material description
  /// @param thickness is the thickness of the material
  MaterialSlab(const Material& material, float thickness);

  /// Scale the material thickness by the given factor.
  /// @param scale Factor by which to scale the thickness
  void scaleThickness(float scale);

  /// Check if the material is vacuum.
  /// @return True if the material is vacuum or thickness is zero/negative
  bool isVacuum() const { return m_material.isVacuum() || m_thickness <= 0; }

  /// Access the (average) material parameters.
  /// @return Reference to the material properties
  constexpr const Material& material() const { return m_material; }
  /// Return the thickness.
  /// @return Material thickness in millimeters
  constexpr float thickness() const { return m_thickness; }
  /// Return the radiation length fraction.
  /// @return Thickness as a fraction of radiation length
  constexpr float thicknessInX0() const { return m_thicknessInX0; }
  /// Return the nuclear interaction length fraction.
  /// @return Thickness as a fraction of nuclear interaction length
  constexpr float thicknessInL0() const { return m_thicknessInL0; }

 private:
  Material m_material;
  float m_thickness = 0.0f;
  float m_thicknessInX0 = 0.0f;
  float m_thicknessInL0 = 0.0f;

  static constexpr auto eps = 2 * std::numeric_limits<float>::epsilon();

  constexpr MaterialSlab(const Material& material, float thickness,
                         [[maybe_unused]] bool dummy)
      : m_material(material), m_thickness(thickness) {
    m_thicknessInX0 = (eps < material.X0()) ? (thickness / material.X0()) : 0;
    m_thicknessInL0 = (eps < material.L0()) ? (thickness / material.L0()) : 0;
  }

  /// @brief Check if two materials are exactly equal.
  ///
  /// This is a strict equality check, i.e. the materials must have identical
  /// properties.
  ///
  /// @param lhs is the left hand side material
  /// @param rhs is the right hand side material
  ///
  /// @return true if the materials are equal
  friend constexpr bool operator==(const MaterialSlab& lhs,
                                   const MaterialSlab& rhs) {
    // t/X0 and t/L0 are dependent variables and need not be checked
    return (lhs.m_material == rhs.m_material) &&
           (lhs.m_thickness == rhs.m_thickness);
  }
};

/// Stream operator for MaterialSlab
/// @param os Output stream
/// @param materialSlab MaterialSlab to output
/// @return Reference to output stream
std::ostream& operator<<(std::ostream& os, const MaterialSlab& materialSlab);

/// @brief Type alias for a vector of material slabs
/// @details Used to store a collection of material slabs in sequence
using MaterialSlabVector = std::vector<MaterialSlab>;

/// @brief Type alias for a matrix of material slabs
/// @details Used to store a 2D collection of material slabs
using MaterialSlabMatrix = std::vector<MaterialSlabVector>;

/// list of point used in the mapping of a volume
using RecordedMaterialVolumePoint =
    std::vector<std::pair<Acts::MaterialSlab, std::vector<Acts::Vector3>>>;

}  // namespace Acts

#if defined(__clang__) && defined(__x86_64__)
#pragma float_control(pop)
#endif

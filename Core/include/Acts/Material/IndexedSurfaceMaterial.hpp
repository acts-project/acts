// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <vector>

namespace Acts {

struct MaterialSlabAccessor {
  using material_type = MaterialSlab;
  // Simple forwarding operator
  /// @brief Access the material slab stored in the vector
  inline static const MaterialSlab& slab(const material_type& mslab) {
    return mslab;
  }
  // Simple forwarding operator
  /// @brief Access the material slab stored in the vector
  inline static MaterialSlab& slab(material_type& mslab) { return mslab; }
};

/// @class IndexedSurfaceMaterial
///
/// It extends the @c ISurfaceMaterial base class and allows to create
/// material maps associated to a grid structure
///
/// @tparam grid_type is the type of the grid used here
/// @tparam material_accessor_type is the type of the accessor to the material
///
/// It is templated on the material type and a slab accessor type in order
/// to allow it to be used in the material recording as well.
template <typename grid_type,
          typename material_accessor_type = MaterialSlabAccessor>
class IndexedSurfaceMaterial : public ISurfaceMaterial {
 public:
  /// @brief Constructor for indexed surface material
  /// @param material the material vector in a flat format
  /// @param grid the index grid steering the access to the material vector
  /// @param gcasts global casts -> casts a grid position from global
  /// @param laccessors local accessors -> accessors to the local grid
  /// @param transform  transform from global frame into map frame
  IndexedSurfaceMaterial(
      std::vector<typename material_accessor_type::material_type>&& material,
      grid_type&& grid, const std::vector<BinningValue>& gcasts,
      const std::vector<std::size_t>& laccessors,
      const Transform3& transform = Transform3::Identity())
      : m_material(std::move(material)),
        m_grid(std::move(grid)),
        m_globalCasts(gcasts),
        m_localAccessors(laccessors),
        m_transform(transform) {
    if (gcasts.size() != grid_type::DIM) {
      throw std::invalid_argument(
          "The number of casts must match the grid dimension");
    }

    if (gcasts.size() != laccessors.size()) {
      throw std::invalid_argument(
          "The number of casts must match the number of local accessors");
    }
  }

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector2&) const
  const MaterialSlab& materialSlab(const Vector2& lp) const final {
    // Access via local position lookup
    std::size_t index = m_grid.atPosition(
        GridAccessHelpers::accessLocal<grid_type>(lp, m_localAccessors));
    return material_accessor_type::slab(m_material[index]);
  }

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector3&) const
  const MaterialSlab& materialSlab(const Vector3& gp) const final {
    // Access via (transformed) global position lookup
    auto index = m_grid.atPosition(GridAccessHelpers::castPosition<grid_type>(
        m_transform * gp, m_globalCasts));
    return material_accessor_type::slab(m_material[index]);
  }

  /// @copydoc ISurfaceMaterial::materialSlab(std::size_t bin0, std::size_t bin1) const
  ///
  /// The ISurface material class expects bins in [ 0 - n ] x [ 0, m ], we will
  /// convert into this format, but it is gonna be slow
  const MaterialSlab& materialSlab(std::size_t bin0,
                                   std::size_t bin1) const final {
    // Access via bin0 and bin1
    Acts::Vector2 lposition =
        GridAccessHelpers::toLocal(m_grid, bin0, bin1, m_localAccessors[0u]);
    return materialSlab(lposition);
  }

  /// Scale operator
  ///
  /// @param scale is the scale factor applied
  ISurfaceMaterial& operator*=(ActsScalar scale) final {
    for (unsigned int im = 0; im < m_material.size(); ++im) {
      material_accessor_type::slab(m_material[im]).scaleThickness(scale);
    }
    return (*this);
  }

  /// Output Method for std::ostream, to be overloaded by child classes
  std::ostream& toStream(std::ostream& sl) const final {
    sl << "IndexedSurfaceMaterial - with an grid<std::size> material store.";
    return sl;
  }

 private:
  /// @brief The stored material
  std::vector<typename material_accessor_type::material_type> m_material;

  /// @brief The grid
  grid_type m_grid;

  /// @brief The casts from global to local
  std::vector<BinningValue> m_globalCasts;

  /// @brief Local accessors - empty assumes direct translation
  std::vector<std::size_t> m_localAccessors;

  /// @brief The transform to the local frame
  Transform3 m_transform;
};

}  // namespace Acts

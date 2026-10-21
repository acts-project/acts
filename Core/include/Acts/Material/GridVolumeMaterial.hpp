// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"

#include <memory>
#include <variant>
#include <vector>

namespace Acts {

/// Concrete volume material backed by a three-dimensional grid.
///
/// Like GridSurfaceMaterial, the binning and storage are independent. Each bin
/// contains constant Material rather than an interpolated parameter vector.
/// Binning is independent of the tracking geometry. Directions on all three
/// axes select coordinates in the material frame, e.g. (r, phi, z); if none are
/// specified, (x, y, z) is used. Partial or duplicate directions are rejected.
/// Open axes have explicit exterior bins, bound axes clamp, and closed axes
/// wrap coordinates periodically. Storage includes under- and overflow bins.
class GridVolumeMaterial final : public IVolumeMaterial {
 public:
  /// One material per bin, including under- and overflow bins
  using Direct = std::vector<Material>;

  /// Per-bin indices into locally owned material
  struct Indexed {
    /// One index per bin, including under- and overflow
    std::vector<std::size_t> indices;
    /// Locally owned material palette
    std::vector<Material> material;
  };

  /// Per-bin indices into shared material
  struct GloballyIndexed {
    /// One index per bin, including under- and overflow
    std::vector<std::size_t> indices;
    /// Shared material palette; keep it unchanged during propagation
    std::shared_ptr<std::vector<Material>> material;
  };

  /// Storage alternatives matching GridSurfaceMaterial
  using Storage = std::variant<Direct, Indexed, GloballyIndexed>;

  /// Construct a grid with fully specified material binning.
  /// @param binning fully specified three-dimensional binning
  /// @param storage material or palette indices for every global bin
  /// @param transform placement of the material frame in global coordinates
  /// @throws std::invalid_argument for inconsistent storage or axis directions
  /// @throws std::domain_error for deferred binning
  GridVolumeMaterial(MultiAxisSpec3D binning, Storage storage,
                     const Transform3& transform = Transform3::Identity());

  /// @copydoc IVolumeMaterial::materialCell
  VolumeMaterialCell materialCell(const Vector3& position) const final;

  /// @copydoc IVolumeMaterial::material(const Vector3&) const
  const Material material(const Vector3& position) const final;

  using IVolumeMaterial::material;

  /// Write the grid binning description.
  /// @param stream destination stream
  /// @return the modified stream
  std::ostream& toStream(std::ostream& stream) const final;

  /// @return material binning specification
  const MultiAxisSpec3D& binning() const { return m_binning; }
  /// @return resolved axes used for material lookup
  const IMultiAxis3D& multiAxis() const { return *m_multiAxis; }
  /// @return material storage backend
  const Storage& storage() const { return m_storage; }
  /// @return placement of the material frame in global coordinates
  const Transform3& transform() const { return m_transform; }

 private:
  const Material& materialAt(std::size_t bin) const;

  MultiAxisSpec3D m_binning;
  std::unique_ptr<IMultiAxis3D> m_multiAxis;
  Storage m_storage;
  Transform3 m_transform;
  std::shared_ptr<const VolumeMaterialCell::Transform> m_globalToGrid;
};

}  // namespace Acts

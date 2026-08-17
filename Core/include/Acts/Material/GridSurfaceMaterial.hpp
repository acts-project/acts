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
#include "Acts/Utilities/IMultiAxis.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"

#include <iosfwd>
#include <memory>
#include <variant>
#include <vector>

namespace Acts {

/// @addtogroup material
/// @{

/// @brief Concrete, non-template surface material backed by a 2D grid.
///
/// The grid geometry is described by a (fully specified) @c MultiAxisSpec2D,
/// resolved internally into an @c IMultiAxis2D that maps a local surface
/// position to a flattened bin index, see @c IMultiAxis. That index then
/// addresses one of three storage backends:
///
/// - @c Direct : a material slab stored directly per bin
/// - @c Indexed : a per-bin index into a locally owned material vector
/// - @c GloballyIndexed : a per-bin index into a (possibly shared) globally
///   owned material vector
///
/// All three backends share the same binning/indexing logic, so the local
/// lookup, scaling and I/O only need to branch on the storage backend, not on
/// the concrete grid or axis types.
///
/// The grid is always 2D: local (bound) lookup is assumed to already be
/// expressed in the grid's own 2D coordinate system, i.e. @c lp[0] maps
/// directly onto axis 0 and @c lp[1] onto axis 1. Only local lookup is
/// supported - there is no global (Vector3) lookup, since that would require
/// a global-to-local coordinate transform.
class GridSurfaceMaterial final : public ISurfaceMaterial {
 public:
  /// Material slabs stored directly, one per bin (including under-/overflow)
  using Direct = std::vector<MaterialSlab>;

  /// Per-bin indices into a locally owned material vector
  struct Indexed {
    /// One index per bin (including under-/overflow), into @c material
    std::vector<std::size_t> indices;
    /// The locally owned material vector, addressed by @c indices
    std::vector<MaterialSlab> material;
  };

  /// Per-bin indices into a (possibly shared) globally owned material vector
  struct GloballyIndexed {
    /// One index per bin (including under-/overflow), into @c material
    std::vector<std::size_t> indices;
    /// The shared material vector, addressed by @c indices
    std::shared_ptr<std::vector<MaterialSlab>> material;
    /// Indicate if entries in @c material are shared with other grids, e.g.
    /// by running a compression/clustering algorithm.
    ///
    /// It is the responsibility of the user to set this flag correctly.
    /// @c scale() throws if it is set, since scaling shared entries would
    /// have an unpredictable outcome.
    bool sharedEntries = false;
  };

  /// The material storage backend: one of @c Direct, @c Indexed or
  /// @c GloballyIndexed
  using Storage = std::variant<Direct, Indexed, GloballyIndexed>;

  /// Construct from a 2D binning spec and a storage backend
  ///
  /// @param binning the fully specified 2D multi-axis binning spec
  /// @param storage the material storage backend
  /// @param splitFactor pre/post splitting directive
  /// @param mappingType type of surface mapping associated to the surface
  /// @throws std::domain_error if @p binning is not fully specified
  /// @throws std::invalid_argument if the storage size does not match the
  ///         number of bins implied by @p binning (including under-/overflow)
  explicit GridSurfaceMaterial(MultiAxisSpec2D binning, Storage storage,
                               double splitFactor = 1.,
                               MappingType mappingType = MappingType::Default);

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector2&) const
  const MaterialSlab& materialSlab(const Vector2& lp) const final;

  /// @copydoc ISurfaceMaterial::materialSlab(const Vector3&) const
  /// @deprecated Global (Vector3) lookup is not supported; use
  ///             materialSlab(const Vector2&) with a prior
  ///             Surface::globalToLocal() call instead.
  [[deprecated(
      "Use materialSlab(const Vector2& lp) with a prior "
      "Surface::globalToLocal() call instead")]] const MaterialSlab&
  materialSlab(const Vector3& gp) const final;

  using ISurfaceMaterial::materialSlab;

  /// @copydoc ISurfaceMaterial::localAxisDirections() const
  std::vector<AxisDirection> localAxisDirections() const final { return {}; }

  /// @copydoc ISurfaceMaterial::scale(double)
  /// @throws std::invalid_argument if the storage is @c GloballyIndexed with
  ///         @c sharedEntries set
  ISurfaceMaterial& scale(double factor) final;

  /// @copydoc ISurfaceMaterial::toStream(std::ostream&) const
  std::ostream& toStream(std::ostream& sl) const final;

  /// Return the 2D multi-axis binning spec
  /// @return const reference to the binning spec
  const MultiAxisSpec2D& binning() const { return m_binning; }

  /// Return the multi-axis used for local position lookup
  /// @return const reference to the resolved multi-axis
  const IMultiAxis2D& multiAxis() const { return *m_multiAxis; }

  /// Return the material storage backend
  /// @return const reference to the storage variant
  const Storage& storage() const { return m_storage; }

 private:
  /// The 2D multi-axis binning spec
  MultiAxisSpec2D m_binning;
  /// The multi-axis built from @c m_binning, doing local position -> index
  std::unique_ptr<IMultiAxis2D> m_multiAxis;
  /// The material storage backend
  Storage m_storage;
};

/// @}

}  // namespace Acts

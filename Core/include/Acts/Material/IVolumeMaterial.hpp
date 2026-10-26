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
#include "Acts/Material/VolumeMaterialCell.hpp"

#include <optional>
#include <utility>

namespace Acts {

class Material;

/// @class IVolumeMaterial
///
/// @ingroup material
///
/// Virtual base class of volume material description
//
/// Material associated with a volume, with constant cells for cached lookup.
class IVolumeMaterial {
 public:
  /// Per-propagation lookup state. A cache may be reused across providers:
  /// changing provider invalidates the stored cell even at the same position.
  struct Cache {
    /// Provider that produced the cached cell
    const IVolumeMaterial* owner = nullptr;
    /// Last retrieved constant material cell
    std::optional<VolumeMaterialCell> cell;
  };

  /// Virtual Destructor
  virtual ~IVolumeMaterial() = default;

  /// Access to actual material
  ///
  /// @param position is the request position for the material call
  /// @return The material properties at the given position
  virtual const Material material(const Vector3& position) const {
    return materialCell(position).material();
  }

  /// Retrieve a constant material cell containing a global position.
  /// @param position global lookup position
  /// @return the cell containing the position
  virtual VolumeMaterialCell materialCell(const Vector3& position) const = 0;

  /// Retrieve material, retaining the cell until the position leaves it.
  /// @param position global lookup position
  /// @param cache per-propagation cache, updated on a cell or provider change
  /// @return material valid until the next update of @p cache
  const Material& material(const Vector3& position, Cache& cache) const {
    if (cache.owner != this || !cache.cell || !cache.cell->isInside(position)) {
      auto cell = materialCell(position);
      cache.cell = std::move(cell);
      cache.owner = this;
    }
    return cache.cell->material();
  }

  /// @brief output stream operator
  ///
  /// Prints information about this object to the output stream using the
  /// virtual IVolumeeMaterial::toStream method
  ///
  /// @return modified output stream object
  friend std::ostream& operator<<(std::ostream& out,
                                  const IVolumeMaterial& vm) {
    vm.toStream(out);
    return out;
  }

  /// Output Method for std::ostream, to be overloaded by child classes
  /// @param sl Output stream to write to
  /// @return Reference to the output stream for method chaining
  virtual std::ostream& toStream(std::ostream& sl) const = 0;
};

}  // namespace Acts

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfacePlacementBase.hpp"

#include <memory>
#include <utility>

namespace Acts {
class Surface;
class PlanarBounds;
class DiscBounds;
class ISurfaceMaterial;
}  // namespace Acts

namespace ActsExamples {

/// @class GenericDetectorElement
///
/// This is a lightweight type of detector element,
/// it simply implements the base class.
///
class GenericDetectorElement : public Acts::SurfacePlacementBase {
 public:
  using identifier_type = unsigned long long;
  using identifier_diff = long long;
  using Identifier = identifier_type;

  /// Broadcast the ContextType
  using ContextType = Acts::GeometryContext;

  /// Constructor for single sided detector element
  /// - bound to a Plane Surface
  ///
  /// @param identifier is the module identifier
  /// @param transform is the transform that element the layer in 3D frame
  /// @param pBounds is the planar bounds for the planar detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  GenericDetectorElement(
      const Identifier identifier, const Acts::Transform3& transform,
      std::shared_ptr<const Acts::PlanarBounds> pBounds, double thickness,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr);

  /// Constructor for single sided detector element
  /// - bound to a Disc Surface
  ///
  /// @param identifier is the module identifier
  /// @param transform is the transform that element the layer in 3D frame
  /// @param dBounds is the planar bounds for the disc like detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  GenericDetectorElement(
      const Identifier identifier, const Acts::Transform3& transform,
      std::shared_ptr<const Acts::DiscBounds> dBounds, double thickness,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr);

  /// Return local to global transform associated with this detector element
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @note this is called from the surface().localToGlobalTransform(gctx) in the PROXY
  /// mode
  const Acts::Transform3& localToGlobalTransform(
      const Acts::GeometryContext& gctx) const override;

  /// Return local to global transform associated with this detector element
  const Acts::Transform3& nominalTransform() const;

  /// Create and return the surface associated with this detector element.
  ///
  /// This method uses @c shared_from_this() so it must only be called after
  /// the element is already managed by a @c std::shared_ptr. The returned
  /// surface takes shared ownership of this detector element.
  ///
  /// @return Shared pointer to the created surface
  std::shared_ptr<Acts::Surface> createSurface() const;

  /// Convenience factory: constructs the element and immediately calls
  /// @c createSurface(), returning both.
  ///
  /// @param identifier Module identifier
  /// @param transform  Placement transform in 3D global space
  /// @param pBounds    Planar bounds of the detector element
  /// @param thickness  Module thickness
  /// @param material   Optional surface material
  /// @return Pair of (shared element, shared surface)
  static std::pair<std::shared_ptr<GenericDetectorElement>,
                   std::shared_ptr<Acts::Surface>>
  create(const Identifier identifier, const Acts::Transform3& transform,
         std::shared_ptr<const Acts::PlanarBounds> pBounds, double thickness,
         std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr);

  /// Convenience factory: constructs the element and immediately calls
  /// @c createSurface(), returning both.
  ///
  /// @param identifier Module identifier
  /// @param transform  Placement transform in 3D global space
  /// @param dBounds    Disc bounds of the detector element
  /// @param thickness  Module thickness
  /// @param material   Optional surface material
  /// @return Pair of (shared element, shared surface)
  static std::pair<std::shared_ptr<GenericDetectorElement>,
                   std::shared_ptr<Acts::Surface>>
  create(const Identifier identifier, const Acts::Transform3& transform,
         std::shared_ptr<const Acts::DiscBounds> dBounds, double thickness,
         std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr);

  /// Return surface associated with this detector element
  /// @deprecated Use @c createSurface() and hold the returned shared_ptr.
  ///             This method will be removed in a future release.
  [[deprecated(
      "Use createSurface() to get a Surface with shared ownership of the "
      "placement; surface() will be removed in a future release")]]
  const Acts::Surface& surface() const override;

  /// Non-cost access to surface associated with this detector element
  /// @deprecated Use @c createSurface() and hold the returned shared_ptr.
  ///             This method will be removed in a future release.
  [[deprecated(
      "Use createSurface() to get a Surface with shared ownership of the "
      "placement; surface() will be removed in a future release")]]
  Acts::Surface& surface() override;

  /// The maximal thickness of the detector element wrt normal axis
  double thickness() const;

  /// The identifier of the detector element
  Identifier identifier() const;
  /// Is the detector element a sensitive element
  bool isSensitive() const override { return true; }

 private:
  // The element identifier
  Identifier m_elementIdentifier;
  /// the transform for positioning in 3D space
  const Acts::Transform3 m_elementTransform;
  /// weak reference back to the surface (owned by external holders via
  /// shared_ptr)
  mutable std::weak_ptr<Acts::Surface> m_elementSurface;
  /// keeps the surface alive when constructed via the deprecated raw-ptr path;
  /// null after createSurface() is called.
  mutable std::shared_ptr<Acts::Surface> m_legacySurface;
  /// the element thickness
  double m_elementThickness;
  /// store either planar or disc bounds
  std::shared_ptr<const Acts::PlanarBounds> m_elementPlanarBounds = nullptr;
  std::shared_ptr<const Acts::DiscBounds> m_elementDiscBounds = nullptr;
  /// optional material, applied when createSurface() is called
  std::shared_ptr<const Acts::ISurfaceMaterial> m_elementMaterial = nullptr;
};

}  // namespace ActsExamples

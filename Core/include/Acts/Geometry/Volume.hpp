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
#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Geometry/VolumePlacementBase.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <iosfwd>
#include <memory>
#include <optional>

namespace Acts {

class VolumeBounds;

/// @class Volume
///
/// It inherits from GeometryObject for geometry identification
///
/// Base class for all volumes inside the tracking realm, it defines the
/// interface for inherited Volume classes regarding the geometrical
/// information.
class Volume : public GeometryObject {
 public:
  /// @brief Type alias for the axis-aligned bounding box of the volume
  /// @details Used to define the spatial extent of the volume in 3D space
  using BoundingBox = AxisAlignedBoundingBox<Volume, double, 3>;

  /// Explicit constructor with shared arguments
  ///
  /// @param transform is the transform to position the volume in 3D space
  /// @param volbounds is the volume boundary definitions
  explicit Volume(const Transform3& transform,
                  std::shared_ptr<VolumeBounds> volbounds) noexcept;
  /// @brief Constructor that connects the volume to an external alignment
  ///        I.e. the volume may move with the alignment of the surfaces
  ///        The placement of the volume is delegated to the positioner
  /// @param positioner: Reference to the object aligning the volume
  /// @param volbounds is the volume boundary definitions
  explicit Volume(VolumePlacementBase& positioner,
                  std::shared_ptr<VolumeBounds> volbounds) noexcept;

  /// @brief Copy Constructor
  /// @param vol is the source volume for the copy
  Volume(const Volume& vol) noexcept;

  /// Copy Constructor with optional shift
  ///
  /// @param vol is the source volume for the copy
  /// @param shift is the optional shift applied as : shift * vol.transform()
  /// @deprecated: Constructor deprecated in favour of shifted(const Transform3& shift) const
  [[deprecated("Use Volume::shifted(const Transform3& shift) const instead.")]]
  Volume(const Volume& vol, const Transform3& shift);

  /// Shift the volume by a transform
  ///
  /// @param shift is the transform to shift the volume by
  /// @return The shifted volume
  Volume shifted(const GeometryContext& gctx, const Transform3& shift) const;

  ~Volume() noexcept override = default;

  /// Assignment operator
  ///
  /// @param vol is the source volume to be copied
  /// @return Reference to this volume for assignment chaining
  Volume& operator=(const Volume& vol) noexcept;

  /// Move assignment operator
  ///
  /// @param other is the other volume to be moved
  /// @return Reference to this volume for assignment chaining
  Volume& operator=(Volume&& other) noexcept = default;

  /// @brief Get the transformation matrix from the local volume frame
  ///        to the global experiment's frame
  /// @param gctx The current geometry context object, e.g. alignment
  /// @return The local to global transformation matrix
  const Transform3& localToGlobalTransform(const GeometryContext& gctx) const;

  /// @brief Get the transformation matrix from the global experiment's
  //         frame to the local volume frame
  /// @param gctx The current geometry context object, e.g. alignment
  /// @return The global to local transformation matrix
  const Transform3& globalToLocalTransform(const GeometryContext& gctx) const;

  /// @brief Get the transform matrix that positions the volume in 3D space
  /// @deprecated: Function deprecated in favour of localToGlobalTransform
  /// @return Const reference to the transform matrix
  [[deprecated(
      "Use localToGlobalTransform(const GeometryContext& gctx) instead.")]]
  const Transform3& transform() const;

  /// @brief Get the inverse transform matrix of the volume
  /// @deprecated: Function deprecated in favour of globalToLocalTransform
  /// @return Const reference to the inverse transform matrix
  [[deprecated(
      "Use globalToLocalTransform(const GeometryContext& gctx) instead.")]]
  const Transform3& itransform() const;

  /// @brief Set the transform matrix for the volume and update internal state
  /// @param transform The new transform matrix to be applied
  void setTransform(const Transform3& transform);

  /// @brief Get the center position of the volume
  /// @param gctx The current geometry context object, e.g. alignment
  /// @return Const reference to the center position vector
  Vector3 center(const GeometryContext& gctx) const;

  /// @brief Get the center position of the volume
  /// @deprecated: Function deprecated in favour of
  ///               center(const GeometryContext& gctx)
  /// @return Const reference to the center position vector
  [[deprecated("Use center(const GeometryContext& gctx) instead.")]]
  const Vector3& center() const;

  /// @brief Get the volume bounds that define the shape of the volume
  /// @return Const reference to the volume bounds object
  const VolumeBounds& volumeBounds() const;

  /// @brief Get mutable access to the volume bounds
  /// @return Reference to the volume bounds object
  VolumeBounds& volumeBounds();

  /// @brief Get shared pointer to the const volume bounds
  /// @return Const shared pointer to the volume bounds object
  std::shared_ptr<const VolumeBounds> volumeBoundsPtr() const;

  /// @brief Get shared pointer to the mutable volume bounds
  /// @return Shared pointer to the volume bounds object
  std::shared_ptr<VolumeBounds> volumeBoundsPtr();

  /// Set volume bounds and update volume bounding boxes implicitly
  /// @param volbounds The volume bounds to be assigned
  void assignVolumeBounds(std::shared_ptr<VolumeBounds> volbounds);

  /// Set the volume bounds and optionally also update the volume transform
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param volbounds The volume bounds to be assigned
  /// @param transform The transform to be assigned, can be optional
  /// @param logger A logger object to log messages
  virtual void update(const GeometryContext& gctx,
                      std::shared_ptr<VolumeBounds> volbounds,
                      std::optional<Transform3> transform = std::nullopt,
                      const Logger& logger = Acts::getDummyLogger());

  /// Construct bounding box for this shape
  /// @param envelope Optional envelope to add / subtract from min/max
  /// @return Constructed bounding box pointing to this volume
  BoundingBox boundingBox(const Vector3& envelope = {0, 0, 0}) const;

  /// Construct oriented bounding box for this shape
  /// @note This will build an oriented bounding box with an
  ///       envelope value of (0.05, 0.05, 0.05)mm
  /// @return Constructed oriented bounding box pointing to this volume
  BoundingBox orientedBoundingBox() const;

  /// Inside() method for checks
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param gpos is the position to be checked
  /// @param tol is the tolerance parameter
  ///
  /// @return boolean indicator if the position is inside
  bool inside(const GeometryContext& gctx, const Vector3& gpos,
              double tol = 0.) const;

  /// Inside() method for checks
  ///
  /// @param gpos is the position to be checked
  /// @param tol is the tolerance parameter
  /// @deprecated: Function deprecated in favour of
  ///               inside(const GeometryContext& gctx, const Vector3& gpos,
  ///               double tol = 0.)
  ///
  /// @return boolean indicator if the position is inside
  [[deprecated(
      "Use inside(const GeometryContext& gctx, const Vector3& gpos, double tol "
      "= 0.) instead.")]]
  bool inside(const Vector3& gpos, double tol = 0.) const;
  /// The binning position method
  /// - as default the center is given, but may be overloaded
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param aDir is the axis direction for the reference position
  /// @return vector 3D that can be used for the binning
  Vector3 referencePosition(const GeometryContext& gctx,
                            AxisDirection aDir) const override;

  /// @brief Compare this volume with another for equality
  /// @param other The other volume to compare with
  /// @return True if the volumes are equal
  bool operator==(const Volume& other) const;

  /// Produces a 3D visualization of this volume
  /// @param helper The visualization helper describing the output format
  /// @param gctx The geometry context
  /// @param viewConfig The view configuration
  void visualize(IVisualization3D& helper, const GeometryContext& gctx,
                 const ViewConfig& viewConfig) const;

  /// @brief Returns the pointer to the associated volume placement (Might be empty)
  VolumePlacementBase* volumePlacement();

  /// @brief Returns the pointer to the associated volume placement (Might be empty)
  const VolumePlacementBase* volumePlacement() const;

  /// @brief Returns whether the volume is alignable (I.e. it's constructed with a volumePlacement)
  bool isAlignable() const;

 private:
  /// @brief Transform matrix that positions the volume in 3D space
  std::unique_ptr<const Transform3> m_transform{};

  /// @brief Inverse of the transform matrix for efficient calculations
  std::unique_ptr<const Transform3> m_itransform{};

  /// @brief Center position of the volume in global coordinates
  Vector3 m_center{Vector3::Zero()};

  /// @brief Volume bounds that define the shape and extent of the volume
  std::shared_ptr<VolumeBounds> m_volumeBounds;
  /// @brief Pointer to the external volume placement that's connected to the alignment
  VolumePlacementBase* m_placement{nullptr};
};

/**Overload of << operator for std::ostream for debug output*/
/// @param sl Output stream
/// @param vol Volume to output
/// @return Reference to output stream
std::ostream& operator<<(std::ostream& sl, const Volume& vol);

}  // namespace Acts

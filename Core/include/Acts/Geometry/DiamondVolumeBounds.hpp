// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/BoundarySurfaceFace.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"

#include <array>
#include <iosfwd>
#include <memory>
#include <ostream>
#include <vector>

namespace Acts {

class RectangleBounds;
class DiamondBounds;

/// @class DiamondVolumeBounds
///
/// Bounds for a polygonal prism shaped Volume, the orientedSurface(...) method
/// creates a vector of 8 surfaces:
/// 2 Diamond Shape Surfaces (see
/// @image html DiamondBounds.svg)
/// and 6 Rectangular Shape Surfaces.
///
///  BoundarySurfaceFace [index]:
///
///  - negativeFaceXY     [0] : Diamond Acts::PlaneSurface,
///                             parallel to \f$ xy \f$ plane at negative \f$ z \f$
///  - positiveFaceXY     [1] : Diamond Acts::PlaneSurface,
///                             parallel to \f$ xy \f$ plane at positive \f$ z \f$
///  - negativeXFaceYZ12  [2] : Rectangular  Acts::PlaneSurface,
///                             parallel to \f$ yz \f$ plane at negative \f$ x \f$ attached at [-x1,-y1] and [-x2,0]
///  - positiveXFaceYZ12  [3] : Rectangular  Acts::PlaneSurface,
///                             parallel to \f$ yz \f$ plane at positive \f$ x \f$ attached at [x1,-y1] and [x2,0]
///  - negativeXFaceYZ23  [4] : Rectangular  Acts::PlaneSurface,
///                             parallel to \f$ yz \f$ plane at negative \f$ x \f$ attached at [-x2,0] and [-x3,y2]
///  - positiveXFaceYZ23  [5] : Rectangular  Acts::PlaneSurface,
///                             parallel to \f$ yz \f$ plane at positive \f$ x \f$ attached at [x2,0] and [x3,y2]
///  - negativeYFaceZX   [6] : Rectangular  Acts::PlaneSurface,
///                             parallel to \f$ zx \f$ plane at negative \f$ y \f$
///  - positiveYFaceZX   [7] : Rectangular  Acts::PlaneSurface,
///                             parallel to \f$ zx \f$ plane at positive \f$ y \f$
///
class DiamondVolumeBounds : public VolumeBounds {
 public:
  /// @enum BoundValues for access / streaming
  enum BoundValues : unsigned int {
    eHalfLengthX1 = 0,  //!< half length in x at positive y1
    eHalfLengthX2 = 1,  //!< half length in x at y=0
    eHalfLengthX3 = 2,  //!< half length in x at negative y2
    eLengthY1 = 3,      //!<  length in positive y1
    eLengthY2 = 4,      //!<  length in negative y2
    eHalfLengthZ = 5,   //!< half length in z
    eAlphaAngle = 6,    //!< angle alpha between negYFaceZX and faceYZ12
    eBetaAngle = 7,     //!< angle beta between FaceYZ12 and FaceYZ23
    eSize               //!< length of the bounds vector
  };

  enum class Face : unsigned int {

    NegativeZFaceXY = 0,
    PositiveZFaceXY = 1,
    NegativeXFaceYZ12 = 2,  // Acts::PlaneSurface attached to [-x1,-y1] and
                            // [-x2,0] at negative x
    PositiveXFaceYZ12 =
        3,  // Acts::PlaneSurface attached to [x1,-y1] and [x2,0] at positive x
    NegativeXFaceYZ23 =
        4,  // Acts::PlaneSurface attached to [-x2,0] and [-x3,y2] at negative x
    PositiveXFaceYZ23 =
        5,  // Acts::PlaneSurface attached to [x2,0] and [x3,y2] at positive x
    NegativeYFaceZX = 6,
    PositiveYFaceZX = 7

  };

  /// Constructor - the polygonal prism boundaries
  ///
  /// @param x1 is the half length in x at negative y1
  /// @param x2 is the half length in x at y=0
  /// @param x3 is the half length in x at positive y2
  /// @param y1 is the positive y extent
  /// @param y2 is the negative y extent
  /// @param halez is the half length in z
  DiamondVolumeBounds(double x1, double x2, double x3, double y1, double y2,
                      double halez) noexcept(false);

  /// Copy constructor
  /// @param other The other DiamondVolumeBounds to copy from
  DiamondVolumeBounds(const DiamondVolumeBounds& other) = default;

  /// Move constructor
  /// @param other The other DiamondVolumeBounds to move from
  DiamondVolumeBounds(DiamondVolumeBounds&& other) = default;

  /// Copy constructor assignment
  /// @param other The other DiamondVolumeBounds to copy from
  DiamondVolumeBounds& operator=(const DiamondVolumeBounds& other) = default;

  /// Move constructor assignment
  /// @param other The other DiamondVolumeBounds to move from
  DiamondVolumeBounds& operator=(DiamondVolumeBounds&& other) = default;

  /// Default destructor
  ~DiamondVolumeBounds() override = default;

  VolumeBounds::BoundsType type() const final {
    return VolumeBounds::BoundsType::eDiamond;
  }

  /// Return the bound values as dynamically sized vector
  ///
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// This method checks if position in the 3D volume frame
  /// is inside the cylinder
  ///
  /// @param pos is the global position to be checked
  /// @param tol is the tolerance applied
  ///
  /// @return boolean indicator if position is inside
  bool inside(const Vector3& pos, double tol = 0.) const final;

  /// Oriented surfaces, i.e. the decomposed boundary surfaces and the
  /// according navigation direction into the volume given the normal
  /// vector on the surface
  ///
  /// @param transform is the 3D transform to be applied to the boundary
  /// surfaces to position them in 3D space
  ///
  /// It will throw an exception if the orientation prescription is not adequate
  ///
  /// @return a vector of surfaces bounding this volume
  std::vector<OrientedSurface> orientedSurfaces(
      const Transform3& transform = Transform3::Identity()) const final;

  /// Construct bounding box for this shape
  /// @param trf Optional transform
  /// @param envelope Optional envelope to add / subtract from min/max
  /// @param entity Entity to associate this bounding box with
  /// @return Constructed bounding box
  Volume::BoundingBox boundingBox(const Transform3* trf = nullptr,
                                  const Vector3& envelope = {0, 0, 0},
                                  const Volume* entity = nullptr) const final;

  /// Output Method for std::ostream
  /// @param os is the output stream
  /// @return Modified ostream for chaining
  std::ostream& toStream(std::ostream& os) const override;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  /// @return The bound value at the specified index
  double get(BoundValues bValue) const { return m_values[bValue]; }

 private:
  /// The bound values stored in an array
  std::array<double, eSize> m_values{};
  /// The face plane surface parallel to YZ at positive X attached to [x1,y1]
  /// and [x2,0]
  std::shared_ptr<Acts::RectangleBounds> m_FaceYZ12Bounds;
  /// The face plane surface parallel to YZ attached to [x2,0] and [x3,y2]
  std::shared_ptr<Acts::RectangleBounds> m_FaceYZ23Bounds;
  // The face plane surface parallel to ZX at negative Y
  std::shared_ptr<Acts::RectangleBounds> m_negYFaceZXBounds;
  // The face plane surface parallel to ZX at positive Y
  std::shared_ptr<Acts::RectangleBounds> m_posYFaceZXBounds;
  /// The face diamond surface parallel to XY
  std::shared_ptr<Acts::DiamondBounds> m_FaceXYBounds;

  // Helper method to build the surface bounds
  void buildSurfaceBounds();

  // Check the input values for consistency,
  // will throw a logic_exception if consistency is not given
  void checkConsistency() noexcept(false);
};

/// @cond
/// @brief Streaming operator for the polygon volume surfaces
inline std::ostream& operator<<(std::ostream& os,
                                DiamondVolumeBounds::Face face) {
  switch (face) {
    case DiamondVolumeBounds::Face::NegativeZFaceXY:
      return os << "NegativeZFaceXY";
    case DiamondVolumeBounds::Face::PositiveZFaceXY:
      return os << "PositiveZFaceXY";
    case DiamondVolumeBounds::Face::NegativeXFaceYZ12:
      return os << "NegativeXFaceYZ12";
    case DiamondVolumeBounds::Face::PositiveXFaceYZ12:
      return os << "PositiveXFaceYZ12";
    case DiamondVolumeBounds::Face::NegativeXFaceYZ23:
      return os << "NegativeXFaceYZ23";
    case DiamondVolumeBounds::Face::PositiveXFaceYZ23:
      return os << "PositiveXFaceYZ23";
    case DiamondVolumeBounds::Face::NegativeYFaceZX:
      return os << "NegativeYFaceZX";
    case DiamondVolumeBounds::Face::PositiveYFaceZX:
      return os << "PositiveYFaceZX";
    default:
      return os << "Unknown";
  }
}
/// @endcond

}  // namespace Acts

// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// Partially abstract base class which provides an interface to visualization
/// helper classes. It provides a number of methods that all the helpers need to
/// conform to. It also provides a `color_type` typedef, but not all of the
/// helpers actually support that.
///
class IVisualization {
 public:
  /// The color typedef. It's an array of three numbers [0, 255] indicating RGB
  /// color values.
  ///
  using color_type = std::array<int, 3>;

  /// The face type
  using face_type = std::vector<size_t>;

  /// Draw a vertex at a given location and a color.
  /// @param vtx The vertex position
  /// @param color The color
  ///
  virtual void vertex(const Vector3D& vtx,
                      color_type color = {120, 120, 120}) = 0;

  /// Draw a face that connects a list of vertices.
  /// @note Depending on the helper implementation, out of plane vertices might
  /// be handled differently.
  /// @param vtxs The vertices that make up the face
  /// @param color The color of the face
  ///
  virtual void face(const std::vector<Vector3D>& vtxs,
                    color_type color = {120, 120, 120}) = 0;

  /// Draw a faces that connects a list of vertices - expert only
  ///
  /// @note Depending on the helper implementation, out of plane vertices might
  /// be handled differently.
  /// @param vtxs The vertices that make up the faceS
  /// @param faces The face presectiotions (i.e. connecting vertices)
  /// @param color The color of the face
  ///
  virtual void faces(const std::vector<Vector3D>& vtxs,
                     const std::vector<face_type>& faces,
                     color_type color = {120, 120, 120}) = 0;

  /// Draw a line from a vertex to another
  /// @param a The start vertex
  /// @param b The end vertex
  /// @param color The color of the line
  ///
  virtual void line(const Vector3D& a, const Vector3D& b,
                    color_type color = {120, 120, 120}) = 0;

  /// Write the content of the helper to an outstream.
  /// @param os The output stream
  ///
  virtual void write(std::ostream& os) const = 0;

  /// Remove all contents of this helper
  ///
  virtual void clear() = 0;

  /// Below are helper functions, which share the same interface as the ones
  /// above, but explicitly accept float values (instead of double), converts
  /// them and calls the above methods.
  ///

  /// @copydoc Acts::IVisualization::vertex(const Vector3D&, color_type)
  ///
  void vertex(const Vector3F& vtx, color_type color = {120, 120, 120}) {
    Vector3D vtxd = vtx.template cast<double>();
    vertex(vtxd, color);
  }

  /// @copydoc Acts::IVisualization::face(std::vector<Vector3F>&, color_type)
  ///
  void face(const std::vector<Vector3F>& vtxs,
            color_type color = {120, 120, 120}) {
    std::vector<Vector3D> vtxsd;
    std::transform(vtxs.begin(), vtxs.end(), std::back_inserter(vtxsd),
                   [](auto& v) { return v.template cast<double>(); });
    face(vtxsd, color);
  }

  ///  @copydoc Acts::IVisualization::line(const Vector3F&, const Vector3F&,
  /// color_type)
  ///
  void line(const Vector3F& a, const Vector3F& b,
            color_type color = {120, 120, 120}) {
    Vector3D ad = a.template cast<double>();
    Vector3D bd = b.template cast<double>();
    line(ad, bd, color);
  }
};

/**
 * Overload of the << operator to facilitate writing to streams.
 * @param os The output stream
 * @param hlp The helper instance
 */
inline std::ostream& operator<<(std::ostream& os, const IVisualization& hlp) {
  hlp.write(os);
  return os;
}
}  // namespace Acts

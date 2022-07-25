// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"

#include <boost/container/small_vector.hpp>

#include <array>
#include <cmath>
#include <limits>

namespace Acts {
class SpacePoint {
  /////////////////////////////////////////////////////////////////////////////////
  // Public methods:
  /////////////////////////////////////////////////////////////////////////////////

 public:
  SpacePoint() = delete;

  /// The only constructor for Acts::SpacePoint. 
  /// Provided position is changed by beamPos.
  /// Provided alignment and variance are added and scaled by standardDeviations.
  /// All locations are in the Seedfinder coordinate system, with the z axis being in the center of the b-field.
  /// @param globalPos 3D point location of the measurement.
  /// @param beamPos x and y locations of the beamSpot w.r.t. the center of the magnetic field
  /// @param variance measurement uncertainty in z and radius.
  /// @param alignment sensor location uncertainty in z and radius.
  /// @param standardDeviations scale factor for variance + alignment. 99.9999% of all measurements fall within 5 standard deviations.
  /// @param sourceLinks vector of Acts::SourceLink* to the measurement(s) corresponding to this SpacePoint
  SpacePoint(const Acts::Vector3& globalPos,
             const Acts::Vector2& beamPos,
             const Acts::Vector2& variance,
             const Acts::Vector2& alignment,
             const float standardDeviations,
             const boost::container::small_vector<const Acts::SourceLink*, 2> sourceLinks);

  SpacePoint(const SpacePoint& sp);

  ~SpacePoint(){
    for (auto* sl : m_sourceLinks){
      delete sl;
    }
  };

  SpacePoint& operator=(
      const SpacePoint&) = delete;

  float x() const { return m_x; }
  float y() const { return m_y; }
  float z() const { return m_z; }
  float radius() const { return m_r; }
  float phi() const { return m_phi; }
  float varianceR() const { return m_varianceR; }
  float varianceZ() const { return m_varianceZ; }
  float deltaR() const { return m_deltaR; }
  float quality() const { return m_quality; }
  float cotTheta() const { return m_cotTheta; }
  void setDeltaR(float deltaR) { m_deltaR = deltaR; }
  void setCotTheta(float cotTheta) { m_cotTheta = cotTheta; }
  void setQuality(float quality) {
    if (quality >= m_quality) {
      m_quality = quality;
    }
  }
  //SpacePoint owns SourceLinks
  const boost::container::small_vector<const Acts::SourceLink*, 2> getSourceLinks() const {return m_sourceLinks;}
  bool operator==(SpacePoint& other);

 protected:
  float m_x;          // x-coordinate   in beam system coordinates
  float m_y;          // y-coordinate   in beam system coordinates
  float m_z;          // z-coordinate   in beam system coordinetes
  float m_r;          // radius         in beam system coordinates
  float m_phi;        // angle around z in beam system coordinates
  float m_varianceR;  // (measurement + alignment uncertainty) * standardDeviations
  float m_varianceZ;  // (measurement + alignment uncertainty) * standardDeviations
  float m_deltaR;     //
  float m_cotTheta = std::numeric_limits<
      float>::quiet_NaN();  // 1/tanTheta estimated from central+this space
                             // point. Its evaluation requires that the space
                             // point is a candidate for triplet search.
  float m_quality = -std::numeric_limits<
      float>::infinity();  // Quality score of the seed the space point is used
                            // for. Quality can be changed if the space point is
                            // used for a better quality seed.
  boost::container::small_vector<const Acts::SourceLink*, 2> m_sourceLinks;
};

inline bool Acts::SpacePoint::operator==(SpacePoint& other){
  const boost::container::small_vector<const Acts::SourceLink*, 2>& otherSL = other.getSourceLinks();

  if (otherSL.size() != m_sourceLinks.size()){
    return false;
  }
  for (size_t i =0; i < m_sourceLinks.size(); i++){
    if (m_sourceLinks[i] != otherSL[i]){
      return false;
    }
  }
  return true;
}

/////////////////////////////////////////////////////////////////////////////////
// Inline methods
/////////////////////////////////////////////////////////////////////////////////

inline SpacePoint::SpacePoint(const Acts::Vector3& globalPos,
                              const Acts::Vector2& beamPos,
                              const Acts::Vector2& variance,
                              const Acts::Vector2& alignment,
                              const float standardDeviations,
                              const boost::container::small_vector<const Acts::SourceLink*, 2> sourceLinks) {
  m_x = globalPos.x() - beamPos.x();
  m_y = globalPos.y() - beamPos.y();
  m_z = globalPos.z();
  m_phi = std::atan2(m_y, m_x);
  m_r = std::sqrt(m_x * m_x + m_y * m_y);
  m_varianceR = (variance.x() + alignment.x()) * standardDeviations;
  m_varianceZ = (variance.y() + alignment.y()) * standardDeviations;
  m_sourceLinks = sourceLinks;
}

/////////////////////////////////////////////////////////////////////////////////
// Copy constructor
/////////////////////////////////////////////////////////////////////////////////

inline SpacePoint::SpacePoint(
    const SpacePoint& sp) = default;

}  // end of namespace Acts

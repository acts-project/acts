// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// RecordedMaterialTrack.hpp, Acts project MaterialPlugins
///////////////////////////////////////////////////////////////////

#pragma once
#include <vector>
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

using RecordedMaterialProperties = std::pair<MaterialProperties, Vector3D>;

/// @class RecordedMaterialTrack
///
/// @brief Holds the material steps along a track, the direction and the
/// starting point of the track.
///
/// The class RecordedMaterialTrack holds a collection of recorded material
/// properties along a track through the detector, the three dimensional
/// global starting position and start momentum.
///
///
/// @note the RecordedMaterialProperties of the track are assumed
/// to be ordered from the starting position along the starting direction
class RecordedMaterialTrack
{
public:
  /// Default constructor
  RecordedMaterialTrack() = default;

  /// Constructor
  ///
  /// @param startPos three dimensional start position of the track (moved)
  /// @param startDir three dimensional start direction of the track (moved)
  /// @param rmps the collection material steps along the track (moved)
  RecordedMaterialTrack(Vector3D                                startPos,
                        Vector3D                                startDir,
                        std::vector<RecordedMaterialProperties> rmps);

  /// Copy constructor
  /// @param mtrecord The source material track
  RecordedMaterialTrack(const RecordedMaterialTrack& mtrecord) = default;

  /// Copy Move constructor
  /// @param mtrecord The source material track
  RecordedMaterialTrack(RecordedMaterialTrack&& mtrecord) = default;

  /// Assignment operator
  /// @param mtrecord The source material track
  RecordedMaterialTrack&
  operator=(const RecordedMaterialTrack& mtrecord)
      = default;

  /// Assignment Move operator
  /// @param mtrecord The source material track
  RecordedMaterialTrack&
  operator=(RecordedMaterialTrack&& mtrecord)
      = default;

  /// Default destructor
  ~RecordedMaterialTrack() = default;

  /// Return method for the position of the track
  const Vector3D&
  position() const;

  /// Return method for the direction of the track
  const Vector3D&
  direction() const;

  ///@return Reference to the accumalted material steps along the track
  const std::vector<RecordedMaterialProperties>&
  recordedMaterialProperties() const;

private:
  /// start position of the material track
  Vector3D m_startPosition;

  /// start direction of the material track
  Vector3D m_startDirection;

  /// the collected material steps along the track
  std::vector<RecordedMaterialProperties> m_recordedMaterialProperties;
};

inline RecordedMaterialTrack::RecordedMaterialTrack(
    Vector3D                                startPos,
    Vector3D                                startDir,
    std::vector<RecordedMaterialProperties> rmps)
  : m_startPosition(std::move(startPos))
  , m_startDirection(std::move(startDir))
  , m_recordedMaterialProperties(std::move(rmps))
{
}

inline const Vector3D&
RecordedMaterialTrack::position() const
{
  return m_startPosition;
}

inline const Vector3D&
RecordedMaterialTrack::direction() const
{
  return m_startDirection;
}

inline const std::vector<RecordedMaterialProperties>&
RecordedMaterialTrack::recordedMaterialProperties() const
{
  return m_recordedMaterialProperties;
}

}  // namespace

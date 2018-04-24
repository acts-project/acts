// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialTrack.hpp, ACTS project MaterialPlugins
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIALPLUGINS_MATERIALTRACK_H
#define ACTS_MATERIALPLUGINS_MATERIALTRACK_H

#include <vector>
#include "Acts/Plugins/MaterialPlugins/MaterialStep.hpp"

namespace Acts {

/// @class MaterialTrack
///
/// @brief Holds the material steps along a track, the direction and the
/// starting point of the track.
///
/// The class MaterialTrack holds a collection of material steps (containg
/// the material of one step and the  position) along a track through the
/// detector, the three dimensional global starting point of the track and the
/// direction expressed in pseudorapidity eta and azimuthal angle phi.
///
/// This class should be used to create material maps consisting of a collection
/// of MaterialTracks through the detector.
///
class MaterialTrack
{
public:
  /// default constructor
  MaterialTrack() = default;

  /// Constructor handing over all the relevant parameters for the material
  /// track
  /// @param startPos three dimensional global start position of the track
  /// @param theta polar angle indicating of the particle direction
  /// @param phi azimuthal angle indicating of the particle direction
  /// @param materialSteps the collection material steps along the track
  /// @param totX0 is the optional total budget in X0
  /// @param totL0 is the optional total budget in L0
  MaterialTrack(const MaterialStep::Position& startPos,
                double                        theta,
                double                        phi,
                std::vector<MaterialStep>     materialSteps,
                double                        totX0 = 0.,
                double                        totL0 = 0.);

  /// Copy constructor
  MaterialTrack(const MaterialTrack& mtrecord);

  /// Default destructor
  ~MaterialTrack() = default;

  /// Assignment operator
  ///
  /// @param mtrechord is the source object
  MaterialTrack&
  operator=(const MaterialTrack& mtrecord);

  ///@return the polar angle theta for the track
  double
  theta() const;

  ///@return returns the azimuthal angle phi for the track
  double
  phi() const;

  ///@return the polar angle theta for the track
  double
  thicknessInX0() const;

  ///@return returns the azimuthal angle phi for the track
  double
  thicknessInL0() const;

#if !defined(__CLING__)
  const Vector3D
  position() const;
#else
  /// return method for the position of the step
  const MaterialStep::Position
  position() const;
#endif

  ///@return returns the accumalted material steps along the track
  std::vector<MaterialStep>
  materialSteps() const;

private:
  /// start position of the material track
  MaterialStep::Position m_startPosition;
  /// pseudorapidity indicating the first coordinate of the direction of the
  /// material track
  double m_theta = 0.;
  /// azimuthal angle phi indicating the second coordinate of the direction of
  /// the material track
  double m_phi = 0.;
  /// total mapped material in X0 (fast access)
  double m_tX0 = 0.;
  /// total mapped material in L0 (fast access)
  double m_tL0 = 0.;
  /// the collected material steps along the track
  std::vector<MaterialStep> m_materialSteps;
};

}  // end of namespace

inline double
Acts::MaterialTrack::theta() const
{
  return m_theta;
}

inline double
Acts::MaterialTrack::phi() const
{
  return m_phi;
}

inline double
Acts::MaterialTrack::thicknessInX0() const
{
  return m_tX0;
}

inline double
Acts::MaterialTrack::thicknessInL0() const
{
  return m_tL0;
}

#if !defined(__CLING__)
inline const Acts::Vector3D
Acts::MaterialTrack::position() const
{
  return Acts::Vector3D(
      m_startPosition.x, m_startPosition.y, m_startPosition.z);
}
#else
inline const Acts::MaterialStep::Position
Acts::MaterialTrack::position() const
{
  return m_startPosition;
}
#endif

inline std::vector<Acts::MaterialStep>
Acts::MaterialTrack::materialSteps() const
{
  return m_materialSteps;
}

#endif  // ACTS_MATERIALPLUGINS_MATERIALTRACK_H

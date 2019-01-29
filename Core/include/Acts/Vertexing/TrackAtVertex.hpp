// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"

namespace Acts {

template<typename InputTrack>
class TrackAtVertex
{
public:
  /// Deleted default constructor
  TrackAtVertex() = delete;

  /// Parameterized constructor
  /// @param chi2perTrack Chi2 of track
  /// @param paramsAtVertex Fitted perigee parameter
  /// @param originalParams Original perigee parameter
  TrackAtVertex(const double&                chi2perTrack,
                const Acts::BoundParameters& fittedParams,
                const InputTrack& originalTrack);

  /// Returns chi2 of track
  double
  chi2() const;

  /// Returns fitted perigee
  const Acts::BoundParameters&
  fittedPerigee() const;

  /// Returns original perigee
  const InputTrack&
  originalTrack() const;

private:
  /// Chi2 of track
  double m_chi2Track;

  /// Fitted perigee
  Acts::BoundParameters m_fittedParams;

  /// Original perigee
  InputTrack m_originalTrack;
};

}  // namespace Acts

#include "Acts/Vertexing/TrackAtVertex.ipp"
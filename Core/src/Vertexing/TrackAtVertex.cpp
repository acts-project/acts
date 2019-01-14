// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/TrackAtVertex.hpp"

Acts::TrackAtVertex::TrackAtVertex(const double&                chi2perTrack,
                                   const Acts::BoundParameters& paramsAtVertex,
                                   const Acts::BoundParameters& originalParams)
  : m_chi2Track(chi2perTrack)
  , m_paramsAtVertex(paramsAtVertex)
  , m_originalParams(originalParams)
{
}

double
Acts::TrackAtVertex::chi2() const
{
  return m_chi2Track;
}

const Acts::BoundParameters&
Acts::TrackAtVertex::fittedPerigee() const
{
  return m_paramsAtVertex;
}

const Acts::BoundParameters&
Acts::TrackAtVertex::originalPerigee() const
{
  return m_originalParams;
}

// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename InputTrack>
Acts::TrackAtVertex<InputTrack>::TrackAtVertex(
    const double&                chi2perTrack,
    const Acts::BoundParameters& paramsAtVertex,
    const InputTrack&            originalTrack)
  : m_chi2Track(chi2perTrack)
  , m_fittedParams(paramsAtVertex)
  , m_originalTrack(originalTrack)
{
}

template <typename InputTrack>
double
Acts::TrackAtVertex<InputTrack>::chi2() const
{
  return m_chi2Track;
}

template <typename InputTrack>
const Acts::BoundParameters&
Acts::TrackAtVertex<InputTrack>::fittedPerigee() const
{
  return m_fittedParams;
}

template <typename InputTrack>
const InputTrack&
Acts::TrackAtVertex<InputTrack>::originalTrack() const
{
  return m_originalTrack;
}

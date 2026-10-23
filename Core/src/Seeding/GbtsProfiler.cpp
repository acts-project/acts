// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsProfiler.hpp"

#include <algorithm>
#include <utility>

namespace Acts::Experimental {

void GbtsProfiler::add(Event&& event) {
  const std::scoped_lock lock(m_mutex);
  m_events.push_back(std::move(event));
}

std::optional<GbtsProfiler::Event> GbtsProfiler::take(
    const std::uint64_t eventNumber) {
  const std::scoped_lock lock(m_mutex);
  const auto it = std::ranges::find(m_events, eventNumber, &Event::eventNumber);
  if (it == m_events.end()) {
    return std::nullopt;
  }
  Event event = std::move(*it);
  m_events.erase(it);
  return event;
}

std::vector<GbtsProfiler::Event> GbtsProfiler::takeAll() {
  const std::scoped_lock lock(m_mutex);
  return std::exchange(m_events, {});
}

std::string_view GbtsProfiler::name(const NodeCut cut) {
  switch (cut) {
    case NodeCut::EndcapClusterWidth:
      return "EndcapClusterWidth";
  }
  return "Unknown";
}

std::string_view GbtsProfiler::name(const DoubletCut cut) {
  switch (cut) {
    case DoubletCut::IsolatedOuterNode:
      return "IsolatedOuterNode";
    case DoubletCut::MinDeltaRadius:
      return "MinDeltaRadius";
    case DoubletCut::StripCrossing:
      return "StripCrossing";
    case DoubletCut::MaxAbsTau:
      return "MaxAbsTau";
    case DoubletCut::InnerTauWindow:
      return "InnerTauWindow";
    case DoubletCut::OuterTauWindow:
      return "OuterTauWindow";
    case DoubletCut::Z0Histogram:
      return "Z0Histogram";
    case DoubletCut::Z0Range:
      return "Z0Range";
    case DoubletCut::ZAtOuterRadius:
      return "ZAtOuterRadius";
    case DoubletCut::Curvature:
      return "Curvature";
    case DoubletCut::MatchBeforeCreate:
      return "MatchBeforeCreate";
    case DoubletCut::EdgeLimit:
      return "EdgeLimit";
  }
  return "Unknown";
}

std::string_view GbtsProfiler::name(const TripletCut cut) {
  switch (cut) {
    case TripletCut::TauRatio:
      return "TauRatio";
    case TripletCut::DeltaPhi:
      return "DeltaPhi";
    case TripletCut::DeltaCurvature:
      return "DeltaCurvature";
    case TripletCut::ValidationDegenerate:
      return "ValidationDegenerate";
    case TripletCut::ValidationD0:
      return "ValidationD0";
    case TripletCut::ValidationPt:
      return "ValidationPt";
    case TripletCut::ValidationHighPtTauRatio:
      return "ValidationHighPtTauRatio";
    case TripletCut::NeighbourLimit:
      return "NeighbourLimit";
  }
  return "Unknown";
}

std::string_view GbtsProfiler::name(const FilterCut cut) {
  switch (cut) {
    case FilterCut::Chi2X:
      return "Chi2X";
    case FilterCut::Chi2Y:
      return "Chi2Y";
    case FilterCut::Curvature:
      return "Curvature";
    case FilterCut::Z0:
      return "Z0";
  }
  return "Unknown";
}

std::string_view GbtsProfiler::name(const HeadOutcome outcome) {
  switch (outcome) {
    case HeadOutcome::Candidate:
      return "Candidate";
    case HeadOutcome::ChainTooLong:
      return "ChainTooLong";
    case HeadOutcome::ChainTooShort:
      return "ChainTooShort";
    case HeadOutcome::FilterRejected:
      return "FilterRejected";
    case HeadOutcome::NoPathSurvived:
      return "NoPathSurvived";
  }
  return "Unknown";
}

std::string_view GbtsProfiler::name(const CandidateFate fate) {
  switch (fate) {
    case CandidateFate::Seed:
      return "Seed";
    case CandidateFate::Split:
      return "Split";
    case CandidateFate::Clone:
      return "Clone";
    case CandidateFate::TooFewSpacePoints:
      return "TooFewSpacePoints";
  }
  return "Unknown";
}

}  // namespace Acts::Experimental

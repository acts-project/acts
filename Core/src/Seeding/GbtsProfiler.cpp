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

namespace {

/// The names of a label type in the order of its values: the first value
/// without a name is past the last.
template <typename Label>
std::vector<std::string_view> allNames() {
  std::vector<std::string_view> names;
  for (std::uint8_t value = 0;
       GbtsProfiler::name(static_cast<Label>(value)) != "Unknown"; ++value) {
    names.push_back(GbtsProfiler::name(static_cast<Label>(value)));
  }
  return names;
}

}  // namespace

void GbtsProfiler::setEventNumber(const std::uint64_t eventNumber) {
  const std::scoped_lock lock(m_mutex);
  m_eventNumbers[std::this_thread::get_id()] = eventNumber;
}

void GbtsProfiler::add(Event&& event) {
  const std::scoped_lock lock(m_mutex);
  event.eventNumber = m_eventNumbers[std::this_thread::get_id()];
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

std::string_view GbtsProfiler::name(const detail::GbtsNodeCut cut) {
  switch (cut) {
    case detail::GbtsNodeCut::EndcapClusterWidth:
      return "EndcapClusterWidth";
  }
  return "Unknown";
}

std::string_view GbtsProfiler::name(const detail::GbtsDoubletCut cut) {
  switch (cut) {
    case detail::GbtsDoubletCut::IsolatedOuterNode:
      return "IsolatedOuterNode";
    case detail::GbtsDoubletCut::MinDeltaRadius:
      return "MinDeltaRadius";
    case detail::GbtsDoubletCut::StripCrossing:
      return "StripCrossing";
    case detail::GbtsDoubletCut::MaxAbsTau:
      return "MaxAbsTau";
    case detail::GbtsDoubletCut::InnerTauWindow:
      return "InnerTauWindow";
    case detail::GbtsDoubletCut::OuterTauWindow:
      return "OuterTauWindow";
    case detail::GbtsDoubletCut::Z0Histogram:
      return "Z0Histogram";
    case detail::GbtsDoubletCut::Z0Range:
      return "Z0Range";
    case detail::GbtsDoubletCut::ZAtOuterRadius:
      return "ZAtOuterRadius";
    case detail::GbtsDoubletCut::Curvature:
      return "Curvature";
    case detail::GbtsDoubletCut::MatchBeforeCreate:
      return "MatchBeforeCreate";
    case detail::GbtsDoubletCut::EdgeLimit:
      return "EdgeLimit";
  }
  return "Unknown";
}

std::string_view GbtsProfiler::name(const detail::GbtsTripletCut cut) {
  switch (cut) {
    case detail::GbtsTripletCut::TauRatio:
      return "TauRatio";
    case detail::GbtsTripletCut::DeltaPhi:
      return "DeltaPhi";
    case detail::GbtsTripletCut::DeltaCurvature:
      return "DeltaCurvature";
    case detail::GbtsTripletCut::SumRatio:
      return "SumRatio";
    case detail::GbtsTripletCut::ValidationDegenerate:
      return "ValidationDegenerate";
    case detail::GbtsTripletCut::ValidationD0:
      return "ValidationD0";
    case detail::GbtsTripletCut::ValidationPt:
      return "ValidationPt";
    case detail::GbtsTripletCut::ValidationHighPtTauRatio:
      return "ValidationHighPtTauRatio";
    case detail::GbtsTripletCut::NeighbourLimit:
      return "NeighbourLimit";
  }
  return "Unknown";
}

std::string_view GbtsProfiler::name(const detail::GbtsFilterCut cut) {
  switch (cut) {
    case detail::GbtsFilterCut::Chi2X:
      return "Chi2X";
    case detail::GbtsFilterCut::Chi2Y:
      return "Chi2Y";
    case detail::GbtsFilterCut::Curvature:
      return "Curvature";
    case detail::GbtsFilterCut::Z0:
      return "Z0";
    case detail::GbtsFilterCut::PostFitCurvature:
      return "PostFitCurvature";
  }
  return "Unknown";
}

std::string_view GbtsProfiler::name(const detail::GbtsHeadOutcome outcome) {
  switch (outcome) {
    case detail::GbtsHeadOutcome::Candidate:
      return "Candidate";
    case detail::GbtsHeadOutcome::ChainTooLong:
      return "ChainTooLong";
    case detail::GbtsHeadOutcome::ChainTooShort:
      return "ChainTooShort";
    case detail::GbtsHeadOutcome::FilterRejected:
      return "FilterRejected";
    case detail::GbtsHeadOutcome::NoPathSurvived:
      return "NoPathSurvived";
  }
  return "Unknown";
}

std::string_view GbtsProfiler::name(const detail::GbtsCandidateFate fate) {
  switch (fate) {
    case detail::GbtsCandidateFate::Seed:
      return "Seed";
    case detail::GbtsCandidateFate::Split:
      return "Split";
    case detail::GbtsCandidateFate::Clone:
      return "Clone";
    case detail::GbtsCandidateFate::TooFewSpacePoints:
      return "TooFewSpacePoints";
  }
  return "Unknown";
}

std::vector<std::string_view> GbtsProfiler::nodeCutNames() {
  return allNames<detail::GbtsNodeCut>();
}

std::vector<std::string_view> GbtsProfiler::doubletCutNames() {
  return allNames<detail::GbtsDoubletCut>();
}

std::vector<std::string_view> GbtsProfiler::tripletCutNames() {
  return allNames<detail::GbtsTripletCut>();
}

std::vector<std::string_view> GbtsProfiler::filterCutNames() {
  return allNames<detail::GbtsFilterCut>();
}

std::vector<std::string_view> GbtsProfiler::headOutcomeNames() {
  return allNames<detail::GbtsHeadOutcome>();
}

std::vector<std::string_view> GbtsProfiler::candidateFateNames() {
  return allNames<detail::GbtsCandidateFate>();
}

}  // namespace Acts::Experimental

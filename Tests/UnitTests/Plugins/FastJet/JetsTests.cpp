// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsPlugins/FastJet/Jets.hpp"

using namespace Acts;
using namespace ActsPlugins;

class Track {
 public:
  static constexpr float mass = 139.57061 * UnitConstants::MeV;

  Track(float pt, float eta, float phi) {
    Vector3 p3 = Vector3::Zero();
    p3[0] = pt * std::cos(phi);
    p3[1] = pt * std::sin(phi);
    p3[2] = pt * std::sinh(eta);
    float e = std::sqrt(mass * mass + p3.squaredNorm());
    m_fourMom[0] = p3[0];
    m_fourMom[1] = p3[1];
    m_fourMom[2] = p3[2];
    m_fourMom[3] = e;
  }

  Vector4 fourMomentum() const { return m_fourMom; }

 private:
  Vector4 m_fourMom{};
};

bool operator==(Track const& lhs, Track const& rhs) {
  return lhs.fourMomentum() == rhs.fourMomentum();
}

class TrackContainer {
 public:
  using TrackProxy = Track;

  TrackContainer() = default;
  void insert(Track track) { m_vec.push_back(std::move(track)); }
  std::size_t size() { return m_vec.size(); }

  using ConstTrackProxy = const Track&;
  ConstTrackProxy getTrack(std::size_t i) {
    if (i < size()) {
      return m_vec[i];
    }
    throw std::runtime_error("Too few tracks");
  }

 private:
  std::vector<Track> m_vec{};
};

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(FastJetSuite)

BOOST_AUTO_TEST_CASE(TruthParticleOneJet) {
  ActsFatras::Barcode barcode;
  ActsPlugins::FastJet::TruthJet<TrackContainer> jet(
      Acts::Vector4(100, 0, 0, 100), ActsPlugins::FastJet::JetLabel::Unknown);
  jet.setConstituents(std::vector<ActsFatras::Barcode>{barcode});
  BOOST_CHECK_EQUAL(jet.constituents().size(), 1);
  BOOST_CHECK_EQUAL(jet.constituents()[0], barcode);
}

BOOST_AUTO_TEST_CASE(SingleTrackJet) {
  TrackContainer tracks;
  tracks.insert(Track(100, 0, 0));

  ActsPlugins::FastJet::TrackJet<TrackContainer> jet(
      Acts::Vector4(100, 0, 0, 100), ActsPlugins::FastJet::JetLabel::Unknown);
  std::vector<TrackContainer::TrackProxy> constituents{tracks.getTrack(0)};
  jet.setConstituents(constituents);

  BOOST_CHECK_EQUAL(jet.constituents().size(), 1);
  BOOST_CHECK(jet.constituents()[0] == tracks.getTrack(0));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <Acts/Plugins/FastJet/TrackJets.hpp>

class Track {
 public:
  static constexpr float mass = 139.57061 * Acts::UnitConstants::MeV;

  Track(float pt, float eta, float phi) {
    Acts::Vector3 p3 = Acts::Vector3::Zero();
    p3[0] = pt * std::cos(phi);
    p3[1] = pt * std::sin(phi);
    p3[2] = pt * std::sinh(eta);
    float e = std::sqrt(mass * mass + p3.squaredNorm());
    m_fourMom[0] = p3[0];
    m_fourMom[1] = p3[1];
    m_fourMom[2] = p3[2];
    m_fourMom[3] = e;
  }

  Acts::Vector4 fourMomentum() const { return m_fourMom; }

 private:
  Acts::Vector4 m_fourMom{};
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

BOOST_AUTO_TEST_CASE(SingleTrack) {
  TrackContainer tracks;
  tracks.insert(Track(100, 0, 0));

  Acts::FastJet::InputTracks inputTracks(tracks);

  Acts::FastJet::TrackJetSequence jetSeq =
      Acts::FastJet::TrackJetSequence::create(inputTracks.fourMomenta());
  std::vector<fastjet::PseudoJet> jets = jetSeq.jets();

  BOOST_CHECK_EQUAL(jets.size(), 1);
  BOOST_CHECK_EQUAL(jets[0].constituents().size(), 1);
  BOOST_CHECK_EQUAL(jets[0].constituents()[0].user_index(), 0);
  BOOST_CHECK_CLOSE(jets[0].pt(), 100, 1e-3);
  BOOST_CHECK_CLOSE(jets[0].eta(), 0, 1e-3);
  BOOST_CHECK_CLOSE(jets[0].phi(), 0, 1e-3);
  BOOST_CHECK_CLOSE(jets[0].m(), Track::mass, 1);
}

BOOST_AUTO_TEST_CASE(TwoTracksTwoJets) {
  TrackContainer tracks;
  tracks.insert(Track(100, 0, 0.0));
  tracks.insert(Track(100, 0, std::numbers::pi));

  Acts::FastJet::InputTracks inputTracks(tracks);

  Acts::FastJet::TrackJetSequence jetSeq =
      Acts::FastJet::TrackJetSequence::create(inputTracks.fourMomenta());
  std::vector<fastjet::PseudoJet> jets = jetSeq.jets();

  BOOST_CHECK_EQUAL(jets.size(), 2);

  std::vector<Track> trks_0 = inputTracks.tracksInJet(jets[0]);
  BOOST_CHECK_EQUAL(trks_0.size(), 1);
  BOOST_CHECK(trks_0[0] == tracks.getTrack(0) ||
              trks_0[0] == tracks.getTrack(1));

  std::vector<Track> trks_1 = inputTracks.tracksInJet(jets[1]);
  BOOST_CHECK_EQUAL(trks_1.size(), 1);
  BOOST_CHECK(trks_1[0] == tracks.getTrack(0) ||
              trks_1[0] == tracks.getTrack(1));
  BOOST_CHECK(trks_0[0] != trks_1[0]);
}

BOOST_AUTO_TEST_CASE(TwoTracksOneJet) {
  TrackContainer tracks;
  tracks.insert(Track(100, 0, 0.0));
  tracks.insert(Track(100, 0, 0.2));

  Acts::FastJet::InputTracks inputTracks(tracks);

  Acts::FastJet::TrackJetSequence jetSeq =
      Acts::FastJet::TrackJetSequence::create(inputTracks.fourMomenta());
  std::vector<fastjet::PseudoJet> jets = jetSeq.jets();

  BOOST_CHECK_EQUAL(jets.size(), 1);

  std::vector<Track> trks_0 = inputTracks.tracksInJet(jets[0]);
  BOOST_CHECK_EQUAL(trks_0.size(), 2);
  BOOST_CHECK(trks_0[0] == tracks.getTrack(0) ||
              trks_0[0] == tracks.getTrack(1));
  BOOST_CHECK(trks_0[1] == tracks.getTrack(0) ||
              trks_0[1] == tracks.getTrack(1));
  BOOST_CHECK(trks_0[0] != trks_0[1]);
}

BOOST_AUTO_TEST_CASE(TracksInJetCore) {
  TrackContainer tracks;
  tracks.insert(Track(100, 0, 0));
  tracks.insert(Track(10, 0.05, 0));
  tracks.insert(Track(10, -0.05, 0));
  tracks.insert(Track(10, 0.2, 0));
  tracks.insert(Track(10, -0.2, 0));

  Acts::FastJet::InputTracks inputTracks(tracks);

  Acts::FastJet::TrackJetSequence jetSeq =
      Acts::FastJet::TrackJetSequence::create(inputTracks.fourMomenta());
  std::vector<fastjet::PseudoJet> jets = jetSeq.jets();

  BOOST_REQUIRE_EQUAL(jets.size(), 1);

  std::vector<Track> trks = inputTracks.tracksInJet(jets[0], 0.1);
  BOOST_CHECK_EQUAL(trks.size(), 3);

  BOOST_CHECK(std::find(trks.begin(), trks.end(), tracks.getTrack(0)) !=
              trks.end());
  BOOST_CHECK(std::find(trks.begin(), trks.end(), tracks.getTrack(1)) !=
              trks.end());
  BOOST_CHECK(std::find(trks.begin(), trks.end(), tracks.getTrack(2)) !=
              trks.end());
  BOOST_CHECK(std::find(trks.begin(), trks.end(), tracks.getTrack(3)) ==
              trks.end());
  BOOST_CHECK(std::find(trks.begin(), trks.end(), tracks.getTrack(4)) ==
              trks.end());
}

BOOST_AUTO_TEST_CASE(EmptyTrackContainer) {
  Acts::FastJet::TrackJetSequence jetSeq =
      Acts::FastJet::TrackJetSequence::create(
          std::vector<fastjet::PseudoJet>());
  BOOST_CHECK_EQUAL(jetSeq.jets().size(), 0);
}

BOOST_AUTO_TEST_CASE(InvalidCoreRadius) {
  TrackContainer tracks;
  tracks.insert(Track(100, 0, 0));
  Acts::FastJet::InputTracks inputTracks(tracks);
  Acts::FastJet::TrackJetSequence jetSeq =
      Acts::FastJet::TrackJetSequence::create(inputTracks.fourMomenta());
  BOOST_CHECK_THROW(inputTracks.tracksInJet(jetSeq.jets()[0], -1.0),
                    std::invalid_argument);
}

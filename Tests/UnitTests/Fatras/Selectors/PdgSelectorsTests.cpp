// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/PdgParticle.hpp"
#include "ActsFatras/Selectors/PdgSelectors.hpp"
#include "Dataset.hpp"

namespace {
const auto& detector = Dataset::detector;
const auto& electron = Dataset::centralElectron;
const auto& positron = Dataset::centralPositron;
const auto& muon = Dataset::centralMuon;
const auto& antimuon = Dataset::centralAntiMuon;
const auto& pion = Dataset::centralPion;
}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasKinematicCasts)

BOOST_AUTO_TEST_CASE(PdgSelector) {
  ActsFatras::PdgSelector<Acts::PdgParticle::eElectron> selectElectron;

  BOOST_TEST(selectElectron(detector, electron));
  BOOST_TEST(not selectElectron(detector, positron));
  BOOST_TEST(not selectElectron(detector, muon));
  BOOST_TEST(not selectElectron(detector, antimuon));
  BOOST_TEST(not selectElectron(detector, pion));
}

BOOST_AUTO_TEST_CASE(AbsPdgSelector) {
  ActsFatras::AbsPdgSelector<Acts::PdgParticle::eElectron> selectElectronLike;

  BOOST_TEST(selectElectronLike(detector, electron));
  BOOST_TEST(selectElectronLike(detector, positron));
  BOOST_TEST(not selectElectronLike(detector, muon));
  BOOST_TEST(not selectElectronLike(detector, antimuon));
  BOOST_TEST(not selectElectronLike(detector, pion));
}

BOOST_AUTO_TEST_CASE(PdgExcluder) {
  ActsFatras::PdgExcluder<Acts::PdgParticle::eMuon> excludeMuon;

  BOOST_TEST(excludeMuon(detector, electron));
  BOOST_TEST(excludeMuon(detector, positron));
  BOOST_TEST(not excludeMuon(detector, muon));
  BOOST_TEST(excludeMuon(detector, antimuon));
  BOOST_TEST(excludeMuon(detector, pion));
}

BOOST_AUTO_TEST_CASE(AbsPdgExcluder) {
  ActsFatras::AbsPdgExcluder<Acts::PdgParticle::eMuon> excludeMuonLike;

  BOOST_TEST(excludeMuonLike(detector, electron));
  BOOST_TEST(excludeMuonLike(detector, positron));
  BOOST_TEST(not excludeMuonLike(detector, muon));
  BOOST_TEST(not excludeMuonLike(detector, antimuon));
  BOOST_TEST(excludeMuonLike(detector, pion));
}

BOOST_AUTO_TEST_SUITE_END()

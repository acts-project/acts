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
const auto& electron = Dataset::centralElectron;
const auto& positron = Dataset::centralPositron;
const auto& muon = Dataset::centralMuon;
const auto& antimuon = Dataset::centralAntiMuon;
const auto& pion = Dataset::centralPion;
}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasKinematicCasts)

BOOST_AUTO_TEST_CASE(PdgSelector) {
  ActsFatras::PdgSelector<Acts::PdgParticle::eElectron> selectElectron;
  BOOST_TEST(selectElectron(electron));
  BOOST_TEST(not selectElectron(positron));
  BOOST_TEST(not selectElectron(muon));
  BOOST_TEST(not selectElectron(antimuon));
  BOOST_TEST(not selectElectron(pion));
}

BOOST_AUTO_TEST_CASE(AbsPdgSelector) {
  ActsFatras::AbsPdgSelector<Acts::PdgParticle::eElectron> selectElectronLike;
  BOOST_TEST(selectElectronLike(electron));
  BOOST_TEST(selectElectronLike(positron));
  BOOST_TEST(not selectElectronLike(muon));
  BOOST_TEST(not selectElectronLike(antimuon));
  BOOST_TEST(not selectElectronLike(pion));
}

BOOST_AUTO_TEST_CASE(PdgExcluder) {
  ActsFatras::PdgExcluder<Acts::PdgParticle::eMuon> excludeMuon;
  BOOST_TEST(excludeMuon(electron));
  BOOST_TEST(excludeMuon(positron));
  BOOST_TEST(not excludeMuon(muon));
  BOOST_TEST(excludeMuon(antimuon));
  BOOST_TEST(excludeMuon(pion));
}

BOOST_AUTO_TEST_CASE(AbsPdgExcluder) {
  ActsFatras::AbsPdgExcluder<Acts::PdgParticle::eMuon> excludeMuonLike;
  BOOST_TEST(excludeMuonLike(electron));
  BOOST_TEST(excludeMuonLike(positron));
  BOOST_TEST(not excludeMuonLike(muon));
  BOOST_TEST(not excludeMuonLike(antimuon));
  BOOST_TEST(excludeMuonLike(pion));
}

BOOST_AUTO_TEST_SUITE_END()

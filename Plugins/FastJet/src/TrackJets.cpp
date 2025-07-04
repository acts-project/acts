// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/FastJet/TrackJets.hpp"

namespace Acts::FastJet {

TrackJetSequence TrackJetSequence::create(
    std::vector<fastjet::PseudoJet>& tracks,
    const fastjet::JetDefinition& jetDef) {
  fastjet::ClusterSequence cs(tracks, jetDef);
  return TrackJetSequence(cs);
}

std::vector<fastjet::PseudoJet> TrackJetSequence::jets(float ptMin,
                                                       float etaMax) {
  fastjet::Selector sel_eta = fastjet::SelectorAbsEtaMax(etaMax);
  return sel_eta(m_clusterSeq.inclusive_jets(ptMin));
}

}  // namespace Acts::FastJet

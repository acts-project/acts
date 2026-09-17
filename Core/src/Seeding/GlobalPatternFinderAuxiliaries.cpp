// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/detail/GlobalPatternFinderAuxiliaries.hpp"

namespace Acts::Experimental::detail {


template <HitPayload Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
PatternState<Hit_t, Sector_t, Topology_t>::PatternState(
    const OrderedHit& seed,
    const typename Sector_t::Index_t  expSector,          
    const Config* cfg,
    const Acts::Logger* logger)
        : cfg{cfg},
          m_logger{logger},
          lastInsertedHit{seed},
          prevLayerHit{seed},
          lineAnchorHit{seed},
          patTheta{seed->position.theta()},
          expSect{expSector} {
                
        /** Add the new hit */
        hitsPerGroup[Topology_t::groupIndex(*seed)].push_back(seed);

        /** Update the hit counts */
        if (seed->isPrecision()) {
            nPrecisionLayers++;
        } else {
            nTriggerLayers++;
        }
        if (seed.sp()->measuresPhi()) {
            nPhiLayers++;
        }
        updatePatternPhi();
        needLineUpdate = true;
    }






}
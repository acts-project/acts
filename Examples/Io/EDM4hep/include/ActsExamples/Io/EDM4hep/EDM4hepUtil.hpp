// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <functional>

#include "edm4hep/MCParticle.h"
#include "edm4hep/MutableMCParticle.h"
#include "edm4hep/MutableSimTrackerHit.h"
#include "edm4hep/MutableTrackerHit.h"
#include "edm4hep/MutableTrackerHitPlane.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackerHitPlane.h"

namespace ActsExamples {
namespace EDM4hepUtil {

using MapParticleIdFrom =
    std::function<ActsFatras::Barcode(edm4hep::MCParticle particle)>;
using MapParticleIdTo =
    std::function<edm4hep::MCParticle(ActsFatras::Barcode particleId)>;

using MapGeometryIdFrom =
    std::function<Acts::GeometryIdentifier(std::uint64_t cellId)>;
using MapGeometryIdTo =
    std::function<std::uint64_t(Acts::GeometryIdentifier geometryId)>;

ActsFatras::Particle fromParticle(edm4hep::MCParticle from,
                                  MapParticleIdFrom particleMapper);

void toParticle(const ActsFatras::Particle& from,
                edm4hep::MutableMCParticle to);

ActsFatras::Hit fromSimHit(const edm4hep::SimTrackerHit& from,
                           MapParticleIdFrom particleMapper,
                           MapGeometryIdFrom geometryMapper);

void toSimHit(const ActsFatras::Hit& from, edm4hep::MutableSimTrackerHit to,
              MapParticleIdTo particleMapper, MapGeometryIdTo geometryMapper);

void toMeasurement(const ActsExamples::Measurement& from,
                   edm4hep::MutableTrackerHitPlane to,
                   const Cluster* fromCluster,
                   edm4hep::TrackerHitCollection& toClusters);

}  // namespace EDM4hepUtil
}  // namespace ActsExamples

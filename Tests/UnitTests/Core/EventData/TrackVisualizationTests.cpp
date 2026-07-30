#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/AnyTrackProxy.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/ProxyAccessor.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"
#include "Acts/Visualization/TrackVisualization.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Holders.hpp"

#include <memory>
#include <vector>


namespace {

using namespace Acts::UnitLiterals;
using namespace Acts;
using namespace Acts::HashedStringLiteral;

const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
// Helper to create a test track with some data
template <typename track_container_t>
void fillTestTrack(typename track_container_t::TrackProxy track) {
  auto surface = Acts::Surface::makeShared<PerigeeSurface>(Vector3::Zero());

  // Set reference surface
  track.setReferenceSurface(surface);

  // Set parameters
  auto params = track.parameters();
  params[eBoundLoc0] = 1.0;
  params[eBoundLoc1] = 2.0;
  params[eBoundTime] = 3.0;
  params[eBoundPhi] = 0.5;
  params[eBoundTheta] = std::numbers::pi / 4.0;
  params[eBoundQOverP] = 0.1;

  // Set covariance
  auto cov = track.covariance();
  cov.setIdentity();
  cov(eBoundLoc0, eBoundLoc0) = 0.01;
  cov(eBoundQOverP, eBoundQOverP) = 0.0001;

  // Set statistics
  track.chi2() = 12.5f;
  track.nDoF() = 10u;
  track.nMeasurements() = 8u;
  track.nHoles() = 1u;
  track.nOutliers() = 0u;
  track.nSharedHits() = 2u;

  // Set particle hypothesis
  track.setParticleHypothesis(ParticleHypothesis::pion());
    }
}

BOOST_AUTO_TEST_SUITE(TrackVisualization)
{
    TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};

    auto track = tc.makeTrack();
    fillTestTrack<decltype(tc)>(track);

    trackVis(track);

}

BOOST_AUTO_TEST_SUITE_END()
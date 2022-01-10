#pragma once

#include "Acts/Propagator/MultiEigenStepperLoop.hpp"
#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"
#include "Acts/TrackFitting/GaussianSumFitter.hpp"

using PrecompiledGsfPropagator =
    Acts::Propagator<Acts::MultiEigenStepperLoop<>, Acts::Navigator>;
using BaseGsf = Acts::GaussianSumFitter<PrecompiledGsfPropagator>;

struct PrecompiledGsf : private BaseGsf {
  PrecompiledGsf(PrecompiledGsfPropagator&& prop);

  using source_link_it_t =
      std::vector<Acts::Test::TestSourceLink>::const_iterator;
  using start_parameters_t =
      Acts::MultiComponentBoundTrackParameters<Acts::SinglyCharged>;

  Acts::Result<Acts::KalmanFitterResult> fit(
      source_link_it_t begin, source_link_it_t end,
      const start_parameters_t& sParameters,
      const Acts::GsfOptions& options) const;

  Acts::Result<Acts::KalmanFitterResult> fit(
      source_link_it_t begin, source_link_it_t end,
      const Acts::CurvilinearTrackParameters& sparams,
      const Acts::GsfOptions& options) const {
    start_parameters_t multi_pars(sparams.referenceSurface().getSharedPtr(),
                                  sparams.parameters(), sparams.covariance());

    return this->fit(begin, end, multi_pars, options);
  }
};

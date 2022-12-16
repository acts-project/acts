// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <numeric>
#include <optional>
#include <random>
#include <tuple>

using namespace Acts;

#define ACTS_DOES_NOT_COMPILE_BEGIN(x) \
  {                                    \
#if defined(Enable##x)

#define ACTS_DOES_NOT_COMPILE_END() \
#endif                            \
  }                                 \
  }

int main() {
  {
    ACTS_DOES_NOT_COMPILE_BEGIN(SetThroughConstTrackState)
    constexpr TrackStatePropMask kMask = TrackStatePropMask::Predicted;

    // construct trajectory w/ multiple components
    VectorMultiTrajectory t;

    auto i0 = t.addTrackState(kMask);
    // trajectory bifurcates here into multiple hypotheses
    auto i1b = t.addTrackState(kMask, i0);
    auto i2b = t.addTrackState(kMask, i1b);

    // check const-correctness
    const auto& ct = t;
    for (const auto& p : ct.trackStateRange(i2b)) {
      // mutation in this loop doesn't work: does not compile
      p.predicted() = BoundVector::Random();
    }
    ACTS_DOES_NOT_COMPILE_END()
  }

  {
    // make mutable
    VectorMultiTrajectory t;
    auto i0 = t.addTrackState();

    {
      VectorMultiTrajectory::TrackStateProxy tsp = t.getTrackState(i0);
      static_cast<void>(tsp);
      VectorMultiTrajectory::ConstTrackStateProxy ctsp = t.getTrackState(i0);
      static_cast<void>(ctsp);
    }

    ConstVectorMultiTrajectory ct = t;

    ConstVectorMultiTrajectory ctm{std::move(t)};

    {
      ConstVectorMultiTrajectory::ConstTrackStateProxy ctsp =
          ct.getTrackState(i0);

      ACTS_DOES_NOT_COMPILE_BEGIN(SetPredictedThtoughConstTrackState)
      // doesn't compile:
      ctsp.predictedCovariance().setIdentity();
      ACTS_DOES_NOT_COMPILE_END()
    }

    ACTS_DOES_NOT_COMPILE_BEGIN(ConstTrackStateClear)
    // doesn't compile:
    ct.clear();
    ACTS_DOES_NOT_COMPILE_END()
  }
}

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Plugins/ExaTrkX/GNNTrackFitterCPU.hpp"

BOOST_AUTO_TEST_CASE(instantiate) {
  auto traj = std::make_shared<Acts::VectorMultiTrajectory>();
  auto tc = std::make_shared<Acts::VectorTrackContainer>();
  using TrackContainer =
      Acts::TrackContainer<Acts::VectorTrackContainer,
                           Acts::VectorMultiTrajectory, std::shared_ptr>;
  TrackContainer tracks(tc, traj);

  Acts::GNNTrackFitterCPU<TrackContainer> fitter(
      {}, Acts::getDefaultLogger("blub", {}));

  std::vector<std::vector<int>> candidates;
  std::vector<float> features;
  std::vector<Acts::GeometryIdentifier> geoids;
  std::vector<boost::container::static_vector<Acts::SourceLink, 2>> sls;
  
  Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> exts;
  Acts::GNNTrackFitterCPU<TrackContainer>::Options opts{
    Acts::GeometryContext{}, Acts::MagneticFieldContext{}, Acts::CalibrationContext{}, exts};

  fitter(tracks, candidates, features, geoids, sls, opts);
}

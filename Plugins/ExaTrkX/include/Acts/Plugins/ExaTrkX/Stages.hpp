// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <any>
#include <vector>

namespace Acts {

// TODO maybe replace std::any with some kind of variant<unique_ptr<torch>,
// unique_ptr<onnx>>?
// TODO maybe replace input for GraphConstructionBase with some kind of
// boost::multi_array / Eigen::Array

class GraphConstructionBase : public Loggable {
 public:
  GraphConstructionBase(const Acts::Logger &logger) : Loggable(logger) {}
  virtual std::tuple<std::any, std::any> operator()(
      std::vector<float> &inputValues) = 0;
};

class EdgeClassificationBase : public Loggable {
 public:
  EdgeClassificationBase(const Acts::Logger &logger) : Loggable(logger) {}
  virtual std::tuple<std::any, std::any, std::any> operator()(std::any nodes,
                                                    std::any edges) = 0;
};

class TrackBuildingBase : public Loggable {
 public:
  TrackBuildingBase(const Acts::Logger &logger) : Loggable(logger) {}
  virtual std::vector<std::vector<int>> operator()(
      std::any nodes, std::any edges, std::any edgeWeights, std::vector<int> &spacepointIDs) = 0;
};

}  // namespace Acts

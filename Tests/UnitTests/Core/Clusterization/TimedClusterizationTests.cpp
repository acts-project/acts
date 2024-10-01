// This file is part of the ACTS project.
//
// Copyright (C) 2016-2024 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Clusterization/TimedClusterization.hpp"

namespace Acts::Test {
  
  // Define objects
  using Identifier = std::size_t;
  struct Cell {
    Cell(Identifier identifier, int c, int r, double t)
      : id(identifier), column(c), row(r), time(t) {}
    
    Identifier id{};
    int column{0};
    int row{0};
    int label{-1};
    double time{0.};
  };
  
  struct Cluster {
    std::vector<Identifier> ids{};
  };

  using CellCollection = std::vector<Acts::Test::Cell>;
  using ClusterCollection = std::vector<Acts::Test::Cluster>;
    
  BOOST_AUTO_TEST_CASE(TimedGrid_2D_notime) {
  
  }

}

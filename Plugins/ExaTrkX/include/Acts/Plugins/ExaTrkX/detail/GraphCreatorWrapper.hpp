// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <TTree_hits>
#include <graph>
#include <memory>
#include <string>

template <typename T>
class graph_creator;

template <typename T>
class CUDA_graph_creator;

namespace Acts::detail {

class GraphCreatorWrapperBase {
 public:
  virtual ~GraphCreatorWrapperBase() {}
  virtual graph<float> build(TTree_hits<float> &hits, bool print = false) = 0;
};

class GraphCreatorWrapperCpu : public GraphCreatorWrapperBase {
 public:
  GraphCreatorWrapperCpu(const std::string &path);
  ~GraphCreatorWrapperCpu();

  virtual graph<float> build(TTree_hits<float> &hits,
                             bool print = false) override;

 private:
  std::unique_ptr<graph_creator<float>> m_graphCreator;
};

#ifndef ACTS_EXATRKX_CPUONLY
class GraphCreatorWrapperCuda : public GraphCreatorWrapperBase {
 public:
  GraphCreatorWrapperCuda(const std::string &path, int device, int blocks);
  ~GraphCreatorWrapperCuda();

  virtual graph<float> build(TTree_hits<float> &hits,
                             bool print = false) override;

 private:
  std::unique_ptr<CUDA_graph_creator<float>> m_graphCreator;
};
#endif

}  // namespace Acts::detail

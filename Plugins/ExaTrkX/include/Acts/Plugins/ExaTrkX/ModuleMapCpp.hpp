// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/Stages.hpp"

#include <memory>
#include <string>

template<typename T>
class graph_creator;

namespace Acts {
class ModuleMapCpp : GraphConstructionBase {
public:
    struct Config {
        std::string moduleMapPath;
    };

private:
    Config m_cfg;
    std::unique_ptr<graph_creator<float>> m_graphCreator;

public:
    ModuleMapCpp(const Config &cfg);

    std::tuple<std::any, std::any> operator()(std::vector<float> & inputValues, std::size_t numNodes,  const std::vector<uint64_t> &moduleIds, int deviceHint) override;
};

}

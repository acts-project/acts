// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/DefinitionsJsonConverter.hpp"

void Acts::to_json(nlohmann::json& j, const Direction& direction) {
  j = direction.sign();
}

void Acts::from_json(const nlohmann::json& j, Direction& direction) {
  direction = Direction::fromScalarZeroAsPositive(j.get<int>());
}

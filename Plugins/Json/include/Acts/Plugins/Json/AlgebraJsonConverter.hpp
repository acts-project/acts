// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <nlohmann/json.hpp>

// Custom Json encoder/decoders. Maming is mandated by nlohman::json and thus
// can not match our naming guidelines.

void to_json(nlohmann::json& j, const Acts::Transform3& t);

void from_json(const nlohmann::json& j, Acts::Transform3& t);

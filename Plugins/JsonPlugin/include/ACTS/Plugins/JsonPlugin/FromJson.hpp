// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_PLUGINS_JSONPLUGIN_JSONLOADER_H
#define ACTS_PLUGINS_JSONPLUGIN_JSONLOADER_H

#include "ACTS/Utilities/VariantDataFwd.hpp"
#include "lib/json.hpp"

namespace Acts {

variant_data
from_json(const nlohmann::json& node);

}  // namespace Acts

#endif

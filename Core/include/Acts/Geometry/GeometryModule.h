// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

struct ActsGeometryModuleV1 {
  const char* module_abi_tag;
  const char* user_data_type;
  void* (*build)(const void* user_data, const void* logger);
  void (*destroy)(void* handle);
};

extern "C" const ActsGeometryModuleV1* acts_geometry_module_v1(void);

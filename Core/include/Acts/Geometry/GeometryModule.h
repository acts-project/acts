// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

struct ActsGeometryModuleRequestV1 {
  const char* abi_tag;
  void* context;
  void* user_data;
};

struct ActsGeometryModuleV1 {
  const char* module_abi_tag;
  void* (*build)(const ActsGeometryModuleRequestV1* request);
  void (*destroy)(void* handle);
  const char* (*last_error)(void);
};

extern "C" const ActsGeometryModuleV1* acts_geometry_module_v1(void);

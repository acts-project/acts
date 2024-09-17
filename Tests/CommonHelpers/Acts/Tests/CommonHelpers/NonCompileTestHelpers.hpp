// This file is part of the ACTS project.
//
// Copyright (C) 2023-2024 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#define ACTS_DOES_NOT_COMPILE_SUITE_BEGIN(name) int main() {
#define ACTS_DOES_NOT_COMPILE_SUITE_END() }

#define ACTS_DOES_NOT_COMPILE_BEGIN(name) {
#define ACTS_DOES_NOT_COMPILE_END() }

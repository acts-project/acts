// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#if defined(__cpp_concepts)

#define ACTS_REQUIRES(x) requires(x)
#define ACTS_CONCEPT(x) x
#define ACTS_STATIC_CHECK_CONCEPT(check_concept, check_type) \
  static_assert(check_concept<check_type>,                   \
                #check_type " does not fulfill " #check_concept)

#else

#define ACTS_REQUIRES(x)
#define ACTS_CONCEPT(x) typename
#define ACTS_STATIC_CHECK_CONCEPT(concept, type) \
  static_assert(true, "Dummy assertion")

#endif

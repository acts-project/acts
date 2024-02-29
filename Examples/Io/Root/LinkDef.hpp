// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#ifdef __CLING__

#pragma link C++ class vector < vector < int>> + ;
#pragma link C++ class vector < vector < long>> + ;

#pragma link C++ class vector < vector < unsigned int>> + ;
#pragma link C++ class vector < vector < unsigned long>> + ;

#pragma link C++ class vector < vector < std::int32_t>> + ;
#pragma link C++ class vector < vector < std::int64_t>> + ;

#pragma link C++ class vector < vector < std::uint32_t>> + ;
#pragma link C++ class vector < vector < std::uint64_t>> + ;

#pragma link C++ class vector < vector < float>> + ;
#pragma link C++ class vector < vector < double>> + ;

#endif

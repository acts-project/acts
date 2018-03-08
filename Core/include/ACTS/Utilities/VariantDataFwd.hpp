// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_UTILITIES_VARIANTDATA_FWD_H
#define ACTS_UTILITIES_VARIANTDATA_FWD_H 1

#include <boost/variant/recursive_wrapper.hpp>
#include <boost/variant/variant_fwd.hpp>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace Acts {

class variant_map;
class variant_vector;

using variant_data = boost::variant<int,
                                    double,
                                    std::string,
                                    bool,
                                    boost::recursive_wrapper<variant_map>,
                                    boost::recursive_wrapper<variant_vector>>;
}
#endif

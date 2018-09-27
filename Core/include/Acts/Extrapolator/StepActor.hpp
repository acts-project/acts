// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

struct StepActor
{

  //~ struct this_result
  //~ {
    //~ std::vector<MaterialInteraction> materialInteractions;
  //~ };

  //~ using result_type = this_result;

  template <typename propagator_state_t>
  void
  //~ operator()(propagator_state_t& state, result_type& result) const
  operator()(propagator_state_t& state) const
  {
	  
  }
};
}

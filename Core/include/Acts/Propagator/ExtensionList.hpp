// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/extension_list_implementation.hpp"
#include "Acts/Utilities/detail/Extendable.hpp"
#include "Acts/Utilities/detail/MPL/all_of.hpp"
#include "Acts/Utilities/detail/MPL/has_duplicates.hpp"
//~ #include "Acts/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

template <typename... extensions>
struct ExtensionList : private detail::Extendable<extensions...>
{
private:
  static_assert(not detail::has_duplicates_v<extensions...>,
                "same extension type specified several times");

  using detail::Extendable<extensions...>::tuple;

public:

  using detail::Extendable<extensions...>::get;

	template <typename stepper_state_t>
    void
    evaluatek1(const stepper_state_t& state,
               const Vector3D&        bField,
               Vector3D&              k1)
    {    
		using impl = detail::extension_list_impl<extensions...>;
		impl::evaluatek1(tuple(), state, bField, k1);
    }

    template <typename stepper_state_t>
    bool
    evaluatek2(const stepper_state_t& state,
               const double           half_h,
               const Vector3D&        k1,
               const Vector3D&        bField,
               Vector3D&              k2) const
    {
      return true;
    }

    template <typename stepper_state_t>
    bool
    evaluatek3(const stepper_state_t& state,
               const double           half_h,
               const Vector3D&        k2,
               const Vector3D&        bField,
               Vector3D&              k3) const
    {
      return true;
    }

    template <typename stepper_state_t>
    bool
    evaluatek4(const stepper_state_t& state,
               const double           h,
               const Vector3D&        k3,
               const Vector3D&        bField,
               Vector3D&              k4) const
    {
      return true;
    }

	template<typename stepper_state_t>
    bool
    finalizeStep(stepper_state_t&, const double)
    {
		return true;
	}
	
    bool
    evaluateD(const Vector3D& dir,
              const Vector3D& bField1,
              const Vector3D& bField2,
              const Vector3D& bField3,
              const double    h,
              const Vector3D& k1,
              const Vector3D& k2,
              const Vector3D& k3,
              ActsMatrixD<7, 7>& D) const
    {
      return true;
    }
};

}  // namespace Acts

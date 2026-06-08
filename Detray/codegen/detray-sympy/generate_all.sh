#!/bin/bash

python gen_transport_covariance_to_bound_impl.py ../../core/include/detray/propagator/actors/codegen/covariance_transport.hpp
python gen_full_jacobian.py ../../core/include/detray/propagator/actors/codegen/full_jacobian.hpp
python gen_update_rk_transport_jacobian_impl.py ../../core/include/detray/propagator/codegen/update_rk_transport_jacobian.hpp
python gen_transport_jacobian_types.py ../../core/include/detray/propagator/codegen/transport_jacobian.hpp

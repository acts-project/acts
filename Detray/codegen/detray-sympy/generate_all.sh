#!/bin/bash
#
# Regenerate the detray propagator codegen headers using the shared ACTS sympy
# codegen core (../../../codegen). Run this script from anywhere.
#
# The detray-domain modules (detray_sympy.matrices/checks/output) are imported
# directly from this directory (it is on sys.path[0] when running the gen_*.py
# scripts), so only the shared `codegen` package needs to be installed.

set -euo pipefail

cd "$(dirname "$0")"

ACTS_CODEGEN=../../../codegen

uvrun() {
    uv run --python 3.13 --no-project \
        --with-requirements "${ACTS_CODEGEN}/requirements.txt" \
        --with "${ACTS_CODEGEN}" \
        python "$@"
}

uvrun gen_transport_covariance_to_bound_impl.py ../../core/include/detray/propagator/detail/codegen/covariance_transport.hpp
uvrun gen_full_jacobian.py ../../core/include/detray/propagator/detail/codegen/full_jacobian.hpp
uvrun gen_update_rk_transport_jacobian_impl.py ../../core/include/detray/propagator/detail/codegen/update_rk_transport_jacobian.hpp
uvrun gen_transport_jacobian_types.py ../../core/include/detray/propagator/detail/codegen/transport_jacobian.hpp

#include "PrecompiledGsf.hpp"

PrecompiledGsf::PrecompiledGsf(PrecompiledGsfPropagator &&prop) : BaseGsf(std::move(prop)) {}


Acts::Result<Acts::KalmanFitterResult> PrecompiledGsf::fit(
    PrecompiledGsf::source_link_it_t begin,
    PrecompiledGsf::source_link_it_t end,
    const PrecompiledGsf::start_parameters_t& sParameters,
    const Acts::GsfOptions& options) const {
        return BaseGsf::fit(begin, end, sParameters, options);
    }

target_sources(
    ActsCore
    PRIVATE
        EigenStepperError.cpp
        MultiStepperError.cpp
        NavigatorError.cpp
        SympyStepper.cpp
        NavigationTarget.cpp
        PropagatorError.cpp
        StraightLineStepper.cpp
        detail/MaterialEffectsAccumulator.cpp
        detail/PointwiseMaterialInteraction.cpp
        detail/CovarianceEngine.cpp
        detail/JacobianEngine.cpp
        detail/SympyCovarianceEngine.cpp
        detail/SympyJacobianEngine.cpp
)

include(ActsCodegen)

acts_code_generation(
    ADD_TO_TARGET ActsCore
    KEY Core/src/Propagator/codegen/sympy_cov_math.hpp
)

acts_code_generation(
    ADD_TO_TARGET ActsCore
    KEY Core/src/Propagator/codegen/sympy_stepper_math.hpp
)

acts_code_generation(
    ADD_TO_TARGET ActsCore
    KEY Core/src/Propagator/codegen/sympy_jac_math.hpp
)

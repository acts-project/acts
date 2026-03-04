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
    PYTHON detail/generate_sympy_cov.py
    WITH_REQUIREMENTS ${PROJECT_SOURCE_DIR}/codegen/requirements.txt
    WITH ${PROJECT_SOURCE_DIR}/codegen
    ISOLATED
    OUTPUT codegen/sympy_cov_math.hpp
)

acts_code_generation(
    ADD_TO_TARGET ActsCore
    PYTHON generate_sympy_stepper.py
    WITH_REQUIREMENTS ${PROJECT_SOURCE_DIR}/codegen/requirements.txt
    WITH ${PROJECT_SOURCE_DIR}/codegen
    ISOLATED
    OUTPUT codegen/sympy_stepper_math.hpp
)

acts_code_generation(
    ADD_TO_TARGET ActsCore
    PYTHON detail/generate_sympy_jac.py
    WITH_REQUIREMENTS ${PROJECT_SOURCE_DIR}/codegen/requirements.txt
    WITH ${PROJECT_SOURCE_DIR}/codegen
    ISOLATED
    OUTPUT codegen/sympy_jac_math.hpp
)

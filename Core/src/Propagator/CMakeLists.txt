target_sources(
    ActsCore
    PRIVATE
        EigenStepperError.cpp
        MultiStepperError.cpp
        NavigatorError.cpp
        SympyStepper.cpp
        PropagatorError.cpp
        StraightLineStepper.cpp
        detail/PointwiseMaterialInteraction.cpp
        detail/CovarianceEngine.cpp
        detail/JacobianEngine.cpp
        detail/SympyCovarianceEngine.cpp
        detail/SympyJacobianEngine.cpp
)

include(ActsCodegen)

if(ACTS_USE_SYSTEM_LIBS)
    set(_codegen_offline_flags "OFFLINE")
else()
    set(_codegen_offline_flags "")
endif()

acts_code_generation(
      TARGET SympyCovarianceEngine
      PYTHON detail/generate_sympy_cov.py
      WITH_REQUIREMENTS ${CMAKE_SOURCE_DIR}/codegen/requirements.txt
      WITH ${CMAKE_SOURCE_DIR}/codegen
      ISOLATED
      OUTPUT codegen/sympy_cov_math.hpp
      ${_codegen_offline_flags}
)

acts_code_generation(
      TARGET SympyStepperMath
      PYTHON generate_sympy_stepper.py
      WITH_REQUIREMENTS ${CMAKE_SOURCE_DIR}/codegen/requirements.txt
      WITH ${CMAKE_SOURCE_DIR}/codegen
      ISOLATED
      OUTPUT codegen/sympy_stepper_math.hpp
      ${_codegen_offline_flags}
)

acts_code_generation(
      TARGET SympyJacobianEngine
      PYTHON detail/generate_sympy_jac.py
      WITH_REQUIREMENTS ${CMAKE_SOURCE_DIR}/codegen/requirements.txt
      WITH ${CMAKE_SOURCE_DIR}/codegen
      ISOLATED
      OUTPUT codegen/sympy_jac_math.hpp
      ${_codegen_offline_flags}
)

unset(_codegen_offline_flags)

target_link_libraries(
    ActsCore
    PRIVATE SympyCovarianceEngine SympyStepperMath SympyJacobianEngine
)

# message(FATAL_ERROR "stop")

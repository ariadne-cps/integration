add_library(INTEGRATION_SRC OBJECT
    integrator.cpp
    integrator_accessors.cpp
    validated_integration_step.cpp
    integration_step_adapters.cpp
    validated_step_integrator.cpp
    multifunction_integrator.cpp
    selected_inclusion_integrator.cpp
    inclusion_integrator.cpp
    bounder.cpp
    integrator_interface.cpp
    runge_kutta_integrator.cpp
)

target_include_directories(INTEGRATION_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/integration"
)

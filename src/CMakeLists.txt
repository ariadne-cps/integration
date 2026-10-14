add_library(INTEGRATION_SRC OBJECT
    integrator.cpp
    inclusion_integrator.cpp
    bounder.cpp
    runge_kutta_integrator.cpp
)

target_include_directories(INTEGRATION_SRC PRIVATE
    "${PROJECT_SOURCE_DIR}/include/integration"
)

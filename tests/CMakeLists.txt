include(CTest)

set(UNIT_TESTS
    test_bounder
    test_integrator
    test_runge_kutta_integrator
    test_inclusion_integrator
    test_preconditioned_integrator
    test_integrator_configuration
    test_affine_integrator
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-integration integration)
    add_test(NAME ${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    set_property(TEST ${TEST} PROPERTY LABELS "integration")
endforeach()

add_custom_target(tests)
add_dependencies(tests ${UNIT_TESTS})

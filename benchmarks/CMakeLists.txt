set(INTEGRATION_BENCHMARKS
    benchmark_integrators_maximum_step
    benchmark_integrators_single_step
)

# Benchmarking is an interface-only client dependency. The standalone
# Integration library bundles its production object targets, but not
# Logging objects introduced specifically by Benchmarking.
get_target_property(_benchmark_objects logging ARIADNE_PROJECT_OBJECT_TARGETS)
get_target_property(_integration_objects integration ARIADNE_PROJECT_OBJECT_TARGETS)
if(NOT _benchmark_objects OR NOT _integration_objects)
    message(FATAL_ERROR "Missing Integration/Logging object target registration")
endif()
foreach(_object IN LISTS _integration_objects)
    list(REMOVE_ITEM _benchmark_objects "${_object}")
endforeach()

foreach(BENCHMARK IN LISTS INTEGRATION_BENCHMARKS)
    add_executable(${BENCHMARK} ${BENCHMARK}.cpp)
    target_link_libraries(${BENCHMARK} PRIVATE
        ariadne-integration integration ariadne-benchmarking)
    foreach(_object IN LISTS _benchmark_objects)
        target_sources(${BENCHMARK} PRIVATE $<TARGET_OBJECTS:${_object}>)
    endforeach()
endforeach()

add_dependencies(benchmarks ${INTEGRATION_BENCHMARKS})

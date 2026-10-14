set(INTEGRATION_BENCHMARKS
    benchmark_integrators_maximum_step
    benchmark_integrators_single_step
)

foreach(BENCHMARK ${INTEGRATION_BENCHMARKS})
    add_executable(${BENCHMARK} ${BENCHMARK}.cpp)
    target_link_libraries(${BENCHMARK} PRIVATE ariadne-integration integration)
endforeach()

add_dependencies(benchmarks ${INTEGRATION_BENCHMARKS})

include_directories(${Python3_INCLUDE_DIRS})

get_target_property(PYTHON_COMMON_SOURCE_DIR pyariadne-python-common SOURCE_DIR)

add_library(pyariadne-integration-submodule-obj OBJECT
    integration_submodule.cpp
)

target_include_directories(pyariadne-integration-submodule-obj PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${PYTHON_COMMON_SOURCE_DIR}/include"
    $<TARGET_PROPERTY:pyariadne-kernel,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-integration-submodule-obj PRIVATE
    pyariadne-python-common
    pybind11::module
    integration
)
set_target_properties(pyariadne-integration-submodule-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(pyariadne-integration INTERFACE)
target_sources(pyariadne-integration INTERFACE
    $<TARGET_OBJECTS:pyariadne-integration-submodule-obj>
)
target_include_directories(pyariadne-integration INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<TARGET_PROPERTY:pyariadne-kernel,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-integration INTERFACE
    pyariadne-kernel
)

ariadne_add_python_module(
    NAME integration
    MODULE_SOURCE integration_module.cpp
)

find_library(ARIADNE_INTEGRATION_LIBRARY NAMES ariadne-integration)
find_path(ARIADNE_INTEGRATION_INCLUDE_DIR ariadne-integration.hpp
          PATH_SUFFIXES ariadne-integration)

find_package(Threads REQUIRED)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AriadneIntegration DEFAULT_MSG
    ARIADNE_INTEGRATION_LIBRARY
    ARIADNE_INTEGRATION_INCLUDE_DIR
)

if(AriadneIntegration_FOUND)
    get_filename_component(ARIADNE_INTEGRATION_INCLUDE_PARENT_DIR
        "${ARIADNE_INTEGRATION_INCLUDE_DIR}" DIRECTORY)
    if(NOT TARGET AriadneIntegration::ariadne-integration)
        add_library(AriadneIntegration::ariadne-integration UNKNOWN IMPORTED)
        set_target_properties(AriadneIntegration::ariadne-integration PROPERTIES
            IMPORTED_LOCATION "${ARIADNE_INTEGRATION_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES
                "${ARIADNE_INTEGRATION_INCLUDE_PARENT_DIR};${ARIADNE_INTEGRATION_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "Threads::Threads"
        )
    endif()
endif()

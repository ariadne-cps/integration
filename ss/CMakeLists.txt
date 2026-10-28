if(NOT TARGET function)
    ensure_submodule_initialized(function)
    add_subdirectory(function)
endif()

if(NOT TARGET kernel)
    ensure_submodule_initialized(kernel)
    add_subdirectory(kernel)
endif()

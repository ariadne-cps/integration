if(NOT TARGET function)
    ariadne_add_submodule(function)
endif()

if(NOT TARGET logging)
    ariadne_add_submodule(logging)
endif()

# Testing is linked into unit tests only; reuse a target from transitive dependencies.
if(NOT TARGET ariadne-testing)
    ariadne_add_submodule(testing)
endif()

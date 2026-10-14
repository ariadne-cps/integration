find_library(ARIADNE_INTEGRATION_LIBRARY NAMES ariadne-integration)

find_package(Threads REQUIRED)
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(GMP QUIET gmp)
  pkg_check_modules(MPFR QUIET mpfr)
  pkg_check_modules(CAIRO QUIET cairo)
endif()

find_path(GMP_INCLUDE_DIR gmp.h HINTS ${GMP_INCLUDE_DIRS})
find_path(MPFR_INCLUDE_DIR mpfr.h HINTS ${MPFR_INCLUDE_DIRS})
find_path(CAIRO_INCLUDE_DIR cairo.h HINTS ${CAIRO_INCLUDE_DIRS} PATH_SUFFIXES cairo)
find_library(GMP_LIBRARY NAMES gmp libgmp HINTS ${GMP_LIBRARY_DIRS})
find_library(MPFR_LIBRARY NAMES mpfr libmpfr HINTS ${MPFR_LIBRARY_DIRS})
find_library(CAIRO_LIBRARY NAMES cairo libcairo HINTS ${CAIRO_LIBRARY_DIRS})

find_path(ARIADNE_INTEGRATION_INCLUDE_DIR ariadne-integration.hpp PATH_SUFFIXES ariadne-integration)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AriadneIntegration DEFAULT_MSG
  ARIADNE_INTEGRATION_LIBRARY
  ARIADNE_INTEGRATION_INCLUDE_DIR
  GMP_INCLUDE_DIR
  MPFR_INCLUDE_DIR
  CAIRO_INCLUDE_DIR
  GMP_LIBRARY
  MPFR_LIBRARY
  CAIRO_LIBRARY
)

if(AriadneIntegration_FOUND)
  get_filename_component(ARIADNE_INTEGRATION_INCLUDE_PARENT_DIR ${ARIADNE_INTEGRATION_INCLUDE_DIR} DIRECTORY)
  set(ARIADNE_INTEGRATION_INCLUDE_DIRS
    ${ARIADNE_INTEGRATION_INCLUDE_PARENT_DIR}
    ${ARIADNE_INTEGRATION_INCLUDE_DIR}
    ${MPFR_INCLUDE_DIR}
    ${GMP_INCLUDE_DIR}
    ${CAIRO_INCLUDE_DIR}
  )
  if(NOT TARGET AriadneIntegration::ariadne-integration)
    add_library(AriadneIntegration::ariadne-integration UNKNOWN IMPORTED)
    set_target_properties(AriadneIntegration::ariadne-integration PROPERTIES
      IMPORTED_LOCATION "${ARIADNE_INTEGRATION_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${ARIADNE_INTEGRATION_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${CAIRO_LIBRARY};${MPFR_LIBRARY};${GMP_LIBRARY};Threads::Threads"
    )
  endif()
endif()

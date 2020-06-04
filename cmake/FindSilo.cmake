find_package(PkgConfig QUIET)
pkg_check_modules(PC_Silo silo QUIET)

find_path(Silo_INCLUDE_DIR silo.h pmpio.h HINTS ${PC_Silo_INCLUDE_DIRS})
find_library(Silo_LIBRARY NAMES silo siloh5 HINTS ${PC_Silo_LIBRARY_DIRS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Silo DEFAULT_MSG Silo_INCLUDE_DIR Silo_LIBRARY)

mark_as_advanced(Silo_INCLUDE_DIR Silo_LIBRARY)

if(Silo_INCLUDE_DIR AND Silo_LIBRARY)
  add_library(Silo::silo UNKNOWN IMPORTED)
  set_target_properties(Silo::silo PROPERTIES
    IMPORTED_LOCATION ${Silo_LIBRARY}
    INTERFACE_INCLUDE_DIRECTORIES ${Silo_INCLUDE_DIR})
  if(Silo_LINK_FLAGS)
    set_property(TARGET Silo::silo APPEND_STRING PROPERTY INTERFACE_LINK_LIBRARIES " ${Silo_LINK_FLAGS}")
  endif()
endif()

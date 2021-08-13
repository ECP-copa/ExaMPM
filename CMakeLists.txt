# project settings
cmake_minimum_required(VERSION 3.12)

project(ExaMPM LANGUAGES CXX VERSION 0.1.0)

include(GNUInstallDirs)

# find dependencies
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
find_package(Cabana REQUIRED COMPONENTS Cabana::Cajita Cabana::cabanacore)
if( NOT Cabana_ENABLE_MPI )
  message( FATAL_ERROR "Cabana must be compiled with MPI" )
endif()
if( NOT Cabana_ENABLE_CAJITA )
  message( FATAL_ERROR "Cabana must be compiled with Cajita" )
endif()
find_package(Silo REQUIRED)

# optional dependencies (taken from Cabana)
macro(ExaMPM_add_dependency)
  cmake_parse_arguments(EXAMPM_DEPENDENCY "" "PACKAGE;VERSION" "" ${ARGN})
  find_package( ${EXAMPM_DEPENDENCY_PACKAGE} ${EXAMPM_DEPENDENCY_VERSION} QUIET )
  string(TOUPPER "${EXAMPM_DEPENDENCY_PACKAGE}" EXAMPM_DEPENDENCY_OPTION )
  option(
    ExaMPM_REQUIRE_${EXAMPM_DEPENDENCY_OPTION}
    "Require ExaMPM to build with ${EXAMPM_DEPENDENCY_PACKAGE} support" ${EXAMPM_DEPENDENCY_PACKAGE}_FOUND)
  if(ExaMPM_REQUIRE_${EXAMPM_DEPENDENCY_OPTION})
    find_package( ${EXAMPM_DEPENDENCY_PACKAGE} ${EXAMPM_DEPENDENCY_VERSION} REQUIRED )
  endif()
  set(ExaMPM_ENABLE_${EXAMPM_DEPENDENCY_OPTION} ${${EXAMPM_DEPENDENCY_PACKAGE}_FOUND})
endmacro()

ExaMPM_add_dependency( PACKAGE "ALL" VERSION 0.9.2 )

# find Clang Format
find_package( CLANG_FORMAT 10 )

# library
add_subdirectory(src)

# examples
add_subdirectory(examples)

##---------------------------------------------------------------------------##
## Clang Format
##---------------------------------------------------------------------------##
if(CLANG_FORMAT_FOUND)
  file(GLOB_RECURSE FORMAT_SOURCES src/*.cpp src/*.hpp examples/*.cpp examples/*.hpp)
  add_custom_target(format
    COMMAND ${CLANG_FORMAT_EXECUTABLE} -i -style=file ${FORMAT_SOURCES}
    DEPENDS ${FORMAT_SOURCES})
endif()

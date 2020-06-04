# project settings
cmake_minimum_required(VERSION 3.12)

project(ExaMPM LANGUAGES CXX VERSION 0.1.0)

include(GNUInstallDirs)

# find dependencies
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
find_package(MPI REQUIRED)
find_package(Kokkos 3 REQUIRED)
find_package(Cabana REQUIRED COMPONENTS Cabana::Cajita Cabana::cabanacore)
find_package(Silo REQUIRED)

# library
add_subdirectory(src)

# examples
add_subdirectory(examples)

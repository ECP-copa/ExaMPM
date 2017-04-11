SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  FileIO.hh
  Geometry.hh
  Mesh.hh
  Mesh2d.hh
  Particle.hh
  Ring.hh
  Square.hh
  )

APPEND_SET(SOURCES
  FileIO.cc
  Mesh2d.cc
  Ring.cc
  Square.cc
  )

TRIBITS_ADD_LIBRARY(
  exampm_core
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

TRIBITS_ADD_EXECUTABLE(
  exampm
  NOEXESUFFIX
  NOEXEPREFIX
  SOURCES exampm.cc
  INSTALLABLE
  )

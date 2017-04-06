SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Geometry.hh
  Mesh.hh
  Mesh2d.hh
  Particle.hh
  Square.hh
  )

APPEND_SET(SOURCES
  Mesh2d.cc
  Square.cc
  )

TRIBITS_ADD_LIBRARY(
  exampm_core
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

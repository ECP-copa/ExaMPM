SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Box.hh
  FileIO.hh
  Geometry.hh
  Mesh.hh
  NeoHookeanStress.hh
  Particle.hh
  StressModel.hh
  TensorTools.hh
  )

APPEND_SET(SOURCES
  Box.cc
  FileIO.cc
  Mesh.cc
  NeoHookeanStress.cc
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

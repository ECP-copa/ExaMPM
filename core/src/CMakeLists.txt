SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  BoundaryCondition.hh
  Box.hh
  Geometry.hh
  Mesh.hh
  NeoHookeanStress.hh
  NewtonianViscousStress.hh
  Particle.hh
  ProblemManager.hh
  Ring.hh
  Sheet.hh
  Sphere.hh
  StressModel.hh
  TensorTools.hh
  )

APPEND_SET(SOURCES
  BoundaryCondition.cc
  Box.cc
  Geometry.cc
  Mesh.cc
  NewtonianViscousStress.cc
  NeoHookeanStress.cc
  ProblemManager.cc
  Ring.cc
  Sheet.cc
  Sphere.cc
  )

TRIBITS_ADD_LIBRARY(
  exampm_core
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

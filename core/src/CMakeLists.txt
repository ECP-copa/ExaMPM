SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  BoundaryCondition.hh
  Box.hh
  FileIO.hh
  Geometry.hh
  InfinitesimalStrain.hh
  LagrangianFiniteStrain.hh
  MaterialModel.hh
  Mesh.hh
  NeoHookeanStress.hh
  NewtonianViscousStress.hh
  Particle.hh
  ProblemManager.hh
  Sphere.hh
  StrainModel.hh
  StressModel.hh
  TensorTools.hh
  )

APPEND_SET(SOURCES
  BoundaryCondition.cc
  Box.cc
  FileIO.cc
  Geometry.cc
  InfinitesimalStrain.cc
  LagrangianFiniteStrain.cc
  Mesh.cc
  NewtonianViscousStress.cc
  NeoHookeanStress.cc
  ProblemManager.cc
  Sphere.cc
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

set(HEADERS
  ExaMPM_BoundaryConditions.hpp
  ExaMPM_DenseLinearAlgebra.hpp
  ExaMPM_Mesh.hpp
  ExaMPM_ParticleInit.hpp
  ExaMPM_ProblemManager.hpp
  ExaMPM_Solver.hpp
  ExaMPM_TimeIntegrator.hpp
  ExaMPM_TimeStepControl.hpp
  ExaMPM_Types.hpp
  ExaMPM_VelocityInterpolation.hpp
  )

set(SOURCES
  ExaMPM_Mesh.cpp
  )

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  
add_library(exampm ${SOURCES})

target_link_libraries(exampm
  Cabana::Core
  Cabana::Grid
)

target_include_directories(exampm
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  )

install(TARGETS exampm LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

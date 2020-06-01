include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable( DamBreak dam_break.cpp )
target_link_libraries( DamBreak PRIVATE exampm)
target_include_directories( DamBreak PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

add_executable( FreeFall free_fall.cpp )
target_link_libraries( FreeFall PRIVATE exampm)
target_include_directories( FreeFall PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

install(TARGETS DamBreak FreeFall DESTINATION ${CMAKE_INSTALL_BINDIR})

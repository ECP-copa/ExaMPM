add_executable( DamBreak dam_break.cpp )
target_link_libraries( DamBreak PRIVATE exampm)

add_executable( FreeFall free_fall.cpp )
target_link_libraries( FreeFall PRIVATE exampm)

install(TARGETS DamBreak FreeFall DESTINATION ${CMAKE_INSTALL_BINDIR})

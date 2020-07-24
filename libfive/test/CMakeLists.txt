set(SRCS main.cpp
    mesh.cpp
    util/mesh_checks.cpp
    util/shapes.cpp)
set(LIBS five)

if (BUILD_GUILE_BINDINGS AND GUILE_FOUND)
    set(SRCS ${SRCS} guile.cpp)
    set(LIBS ${LIBS} five-guile)
endif(BUILD_GUILE_BINDINGS AND GUILE_FOUND)

add_executable(libfive-test ${SRCS})
target_link_libraries(libfive-test ${LIBS})
target_include_directories(libfive-test PRIVATE .)

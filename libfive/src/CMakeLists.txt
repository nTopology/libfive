# Add a custom command that produces version.cpp, plus
# a dummy output that's not actually produced, in order
# to force version.cmake to always be re-run before the build
ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
           ${CMAKE_CURRENT_BINARY_DIR}/_version.cpp
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/version.cmake)

add_library(five SHARED
    eval/base.cpp
    eval/eval_interval.cpp
    eval/eval_jacobian.cpp
    eval/eval_array.cpp
    eval/eval_deriv_array.cpp
    eval/eval_deriv.cpp
    eval/eval_feature.cpp
    eval/eval_point.cpp
    eval/tape.cpp
    eval/feature.cpp
    render/discrete/heightmap.cpp
    render/discrete/voxels.cpp
    render/brep/xtree.cpp
    render/brep/contours.cpp
    render/brep/mesh.cpp
    render/brep/marching.cpp
    solve/bounds.cpp
    solve/solver.cpp
    tree/cache.cpp
    tree/opcode.cpp
    tree/template.cpp
    tree/tree.cpp
    libfive.cpp

    ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
)

target_include_directories(five SYSTEM PUBLIC
    ${Boost_INCLUDE_DIRS}
    ${PNG_INCLUDE_DIR}
    ${EIGEN_INCLUDE_DIRS}
)
target_include_directories(five PUBLIC
    ../include
)

set(LIBS ${PNG_LIBRARIES})

if(UNIX AND NOT(APPLE))
    set(LIBS ${LIBS} ${CMAKE_THREAD_LIBS_INIT})
endif(UNIX AND NOT(APPLE))

target_link_libraries(five ${LIBS})
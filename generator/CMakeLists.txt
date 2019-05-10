cmake_minimum_required(VERSION 3.12)

include_directories(
        ${Boost_INCLUDE_DIRS}
        ${MPI_INCLUDE_PATH})
link_libraries(
        n-body-lib
        ${Boost_LIBRARIES}
        ${MPI_LIBRARIES})

add_executable(tree-body
        three_body.cpp)

add_executable(cube
        cube.cpp)


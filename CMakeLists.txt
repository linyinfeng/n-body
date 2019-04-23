cmake_minimum_required(VERSION 3.12)

find_package(MPI REQUIRED)
find_package(Boost COMPONENTS
        mpi
        iostreams
        serialization
        unit_test_framework
        program_options
        REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

add_library(n-body-lib
        src/logging.cpp
        src/communication.cpp
        src/config.hpp
        src/data.hpp
        src/random.hpp
        src/tree.hpp
        src/overloaded.hpp
        src/random_body.hpp
        src/communication.hpp
        src/logging.hpp
        src/space.hpp
        src/physical.hpp)

target_include_directories(n-body-lib PRIVATE
        ${Boost_INCLUDE_DIRS}
        ${MPI_INCLUDE_PATH})
target_link_libraries(n-body-lib PRIVATE
        ${Boost_LIBRARIES}
        ${MPI_LIBRARIES})

add_executable(n-body
        main.cpp)

target_include_directories(n-body PRIVATE
        ${Boost_INCLUDE_DIRS}
        ${MPI_INCLUDE_PATH})
target_link_libraries(n-body PRIVATE
        n-body-lib
        ${Boost_LIBRARIES}
        ${MPI_LIBRARIES})

add_executable(n-body-test
        test/test.cpp
        test/tree.cpp)

target_include_directories(n-body-test PRIVATE src
        ${Boost_INCLUDE_DIRS})
target_link_libraries(n-body-test PRIVATE
        n-body-lib
        ${Boost_LIBRARIES})
target_compile_definitions(n-body-test PRIVATE BOOST_TEST_DYN_LINK)

add_test(test_all n-body-test)
cmake_minimum_required(VERSION 3.12)

find_package(MPI REQUIRED)
find_package(Boost COMPONENTS
        mpi
        serialization
        log
        program_options
        REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS -DBOOST_LOG_DYN_LINK)

add_executable(n-body
        src/main.cpp
        src/config.hpp
        src/data.hpp
        src/random.hpp
        src/tree.hpp)

target_include_directories(n-body PRIVATE ${Boost_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})
target_link_libraries(n-body PRIVATE ${Boost_LIBRARIES} ${MPI_LIBRARIES})

cmake_minimum_required(VERSION 3.5)

project(nerdss)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fsanitize=address -fno-omit-frame-pointer")

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to default, '${default_build_type}'.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
endif()

if(${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please use a different directory (e.g. 'build').")
endif()

# Get system information
cmake_host_system_information(RESULT NUM_CPUS QUERY NUMBER_OF_PHYSICAL_CORES)
cmake_host_system_information(RESULT NUM_THREADS QUERY NUMBER_OF_LOGICAL_CORES)

# Print system information
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR} with ${NUM_CPUS} physical cores and ${NUM_THREADS} logical cores.")

include_directories(include)

option(NERDSS_ENABLE_GTEST "Build Google Test based test targets" OFF)

#add_definitions(-D_DEBUG) # comment this line if no need _DEBUG

set(SOURCE_DIRS
    boundary_conditions
    classes
    error
    math
    parser
    reactions
    system_setup
    trajectory_functions
    io
)

foreach(SOURCE_DIR ${SOURCE_DIRS})
    file(GLOB DIR_SOURCES "src/${SOURCE_DIR}/*.cpp")
    list(APPEND SOURCES ${DIR_SOURCES})
endforeach()

add_executable(nerdss ${SOURCES} EXEs/nerdss.cpp)

# Set up external libraries
find_package(GSL REQUIRED)
target_link_libraries(nerdss GSL::gsl GSL::gslcblas)

# Set up header directories
include_directories(include $(GSL_INCLUDE_DIR))

if(NERDSS_ENABLE_GTEST)
    find_package(GTest REQUIRED)
    if(TARGET GTest::gtest_main)
        set(NERDSS_GTEST_MAIN_TARGET GTest::gtest_main)
    elseif(TARGET GTest::Main)
        set(NERDSS_GTEST_MAIN_TARGET GTest::Main)
    else()
        message(FATAL_ERROR "Google Test was found, but no gtest main target is available.")
    endif()

    add_executable(nerdss_gtest_smoke
        tests/gtest/smoke_test.cpp
    )
    target_link_libraries(nerdss_gtest_smoke
        ${NERDSS_GTEST_MAIN_TARGET}
    )

    add_executable(nerdss_gtest_vector_coord
        tests/gtest/vector_coord_test.cpp
        src/classes/class_Coord.cpp
        src/classes/class_Vector.cpp
    )
    target_include_directories(nerdss_gtest_vector_coord PRIVATE include)
    target_link_libraries(nerdss_gtest_vector_coord
        ${NERDSS_GTEST_MAIN_TARGET}
    )

    add_executable(nerdss_gtest_rotation_math
        tests/gtest/rotation_math_test.cpp
        src/classes/class_Coord.cpp
        src/classes/class_Quat.cpp
        src/classes/class_Vector.cpp
        src/math/matrix_functions.cpp
    )
    target_include_directories(nerdss_gtest_rotation_math PRIVATE include)
    target_link_libraries(nerdss_gtest_rotation_math
        ${NERDSS_GTEST_MAIN_TARGET}
    )

    add_executable(nerdss_gtest_math_functions
        tests/gtest/math_functions_test.cpp
        src/math/math_functions.cpp
    )
    target_include_directories(nerdss_gtest_math_functions PRIVATE include)
    target_link_libraries(nerdss_gtest_math_functions
        ${NERDSS_GTEST_MAIN_TARGET}
    )

    add_executable(nerdss_gtest_association_angle_helpers
        tests/gtest/association_angle_helpers_test.cpp
        src/classes/class_Coord.cpp
        src/classes/class_Quat.cpp
        src/classes/class_Vector.cpp
        src/reactions/angleSignIsCorrect.cpp
        src/reactions/areParallel.cpp
        src/reactions/areSameAngle.cpp
        src/reactions/create_arbitrary_vector.cpp
        src/reactions/requiresSignFlip.cpp
    )
    target_include_directories(nerdss_gtest_association_angle_helpers PRIVATE include)
    target_link_libraries(nerdss_gtest_association_angle_helpers
        ${NERDSS_GTEST_MAIN_TARGET}
    )
endif()

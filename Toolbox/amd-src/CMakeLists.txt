cmake_minimum_required(VERSION 3.0)

project(amd Fortran)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

file(GLOB AMD_SRC "*.f")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
    link_libraries(pthread gfortran quadmath openblas)
else() # MSVC COMPILER
    link_libraries(libopenblas)
endif()

if(BUILD_SHARED)
    add_library(${PROJECT_NAME} SHARED ${AMD_SRC})
else()
    add_library(${PROJECT_NAME} STATIC ${AMD_SRC})
endif()

message("AMD Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
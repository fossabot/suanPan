cmake_minimum_required(VERSION 3.0)

project(spmm Fortran)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

file(GLOB SPMM_SRC "*.f")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
    link_libraries(pthread gfortran quadmath openblas)
else() # MSVC COMPILER
    link_libraries(libopenblas)
endif()

if(BUILD_SHARED)
    add_library(${PROJECT_NAME} SHARED ${SPMM_SRC})
else()
    add_library(${PROJECT_NAME} STATIC ${SPMM_SRC})
endif()

message("SPMM Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
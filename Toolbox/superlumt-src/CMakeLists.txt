cmake_minimum_required(VERSION 3.0)

project(superlumt C)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

file(GLOB SUPERLUMT_SRC "*.c")

add_definitions(-DNDEBUG -DPRNTlevel=0 -DAdd_)

if(CMAKE_C_COMPILER_ID MATCHES "GNU") # GNU GCC COMPILER
    add_definitions(-D__PTHREAD)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
    set_source_files_properties(await.c PROPERTIES COMPILE_FLAGS "-O0")
    link_libraries(pthread gfortran quadmath openblas)
elseif(CMAKE_C_COMPILER_ID MATCHES "MSVC")
    add_definitions(-D__OPENMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0 /openmp")
    set_source_files_properties(await.c PROPERTIES COMPILE_FLAGS "/Od")
    link_libraries(libopenblas)
endif()

if(BUILD_SHARED)
    add_library(${PROJECT_NAME} SHARED ${SUPERLUMT_SRC})
else()
    add_library(${PROJECT_NAME} STATIC ${SUPERLUMT_SRC})
endif()

message("SUPERLUMT C_FLAGS: ${CMAKE_C_FLAGS}")
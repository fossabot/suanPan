cmake_minimum_required(VERSION 3.0)

project(superlu C)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

file(GLOB SUPERLU_SRC "*.c")

add_definitions(-DNDEBUG -DPRNTlevel=0 -DAdd_)

if(CMAKE_C_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
    link_libraries(pthread gfortran quadmath openblas)
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0")
    link_libraries(libopenblas)
endif()

if(BUILD_SHARED)
    add_library(${PROJECT_NAME} SHARED ${SUPERLU_SRC})
else()
    add_library(${PROJECT_NAME} STATIC ${SUPERLU_SRC})
endif()

message("SUPERLU C_FLAGS: ${CMAKE_C_FLAGS}")
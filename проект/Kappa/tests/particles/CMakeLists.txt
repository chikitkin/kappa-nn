cmake_minimum_required(VERSION 2.8)
project(particleTest)

# ------------------------------------------------------------------------------
add_definitions(-std=c++11)

# ------------------------------------------------------------------------------
include_directories(    ${KappaSrc}/include
                        ${KappaSrc}
                        )
link_directories(       ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
                        )
set(externLIB       kappa
                    yaml-cpp
                    )
set(SRC  
                particleTest.cpp
)
add_library(testsObjects OBJECT ${SRC})

add_executable(particleTest  $<TARGET_OBJECTS:testsObjects>)
target_link_libraries(particleTest  ${externLIB})















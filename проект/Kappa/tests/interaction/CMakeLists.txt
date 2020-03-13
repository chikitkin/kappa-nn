cmake_minimum_required(VERSION 2.8)
project(interactionTest)

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
                interactionTest.cpp
)
add_library(interactionTestObjects OBJECT ${SRC})

add_executable(interactionTest  $<TARGET_OBJECTS:interactionTestObjects>)
target_link_libraries(interactionTest  ${externLIB})















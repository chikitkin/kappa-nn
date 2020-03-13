cmake_minimum_required(VERSION 2.8)
project(csv2YamlParticles)

# ------------------------------------------------------------------------------
add_definitions(-std=c++11)

# ------------------------------------------------------------------------------
include_directories(    ${KappaSrc}/include
                        )
link_directories(       ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
                        )
set(externLIB       kappa
                    yaml-cpp
                    )
set(SRC  
                csv2YamlParticles.cpp
)
add_library(ObjectsParticles OBJECT ${SRC})

add_executable(csv2YamlParticles  $<TARGET_OBJECTS:ObjectsParticles>)
target_link_libraries(csv2YamlParticles  ${externLIB})















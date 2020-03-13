cmake_minimum_required(VERSION 2.8)
project(kappaLib)

# ------------------------------------------------------------------------------
add_definitions(-std=c++11 -fPIC)

# ------------------------------------------------------------------------------
include_directories(    ${KappaSrc}/include
                        ${KappaSrc}
                        )

set(SRC  
                Particles/particle.cpp
                Particles/atom.cpp
                Particles/molecule.cpp
                Interaction/interaction.cpp
)

add_library(kappaObjects OBJECT ${SRC})

add_library(kappa SHARED $<TARGET_OBJECTS:kappaObjects> )
#add_library(kappa STATIC $<TARGET_OBJECTS:kappaObjects> )











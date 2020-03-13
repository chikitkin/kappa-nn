cmake_minimum_required(VERSION 2.8)
project(kappa)


# Global
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(KappaSrc ${CMAKE_SOURCE_DIR}/src)



add_subdirectory(src)
add_subdirectory(tests/particles)
add_subdirectory(tests/interaction)
add_subdirectory(utils/convert/csv2YamlParticles)













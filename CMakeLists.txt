cmake_minimum_required(VERSION 3.16)

project(SW3_block_and_screen LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)

add_library(sw3_lib
    src/BoundingBox.cpp
    src/ParticleBlock.cpp
    src/BVH8.cpp
    src/StillingerWeber.cpp
)

target_include_directories(sw3_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_options(sw3_lib PRIVATE -O3)

add_executable(sw3_bench
    src/main.cpp
)

target_link_libraries(sw3_bench PRIVATE sw3_lib)

target_include_directories(sw3_bench
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Optional: enable warnings for development builds.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(sw3_bench PRIVATE -Wall -Wextra -Wpedantic -O2 )
endif()

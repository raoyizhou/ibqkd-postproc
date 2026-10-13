cmake_minimum_required(VERSION 3.16)
project(libqkd-postproc VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(MSVC)
  add_compile_options(/W4 /O2)
else()
  add_compile_options(-Wall -Wextra -O3 -mpopcnt)
endif()

add_library(qkdpp STATIC
  src/bitvec.cpp
  src/toeplitz.cpp
  src/cascade.cpp
  src/gf2poly.cpp
)
target_include_directories(qkdpp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# GF(2) 多项式乘法的硬件加速：x86 PCLMULQDQ
# 注意：-mpclmul 只决定「可否发射该指令」，是否真正启用由 src/gf2poly.cpp 的
# 运行时 CPUID 探测决定（两个条件同时成立才走硬件路径），因此在不支持
# PCLMULQDQ 的 CPU 上运行也不会 SIGILL，而是自动回落软件实现。
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mpclmul" QKDPP_HAS_PCLMUL)
if(QKDPP_HAS_PCLMUL)
  target_compile_options(qkdpp PUBLIC -mpclmul)
  message(STATUS "PCLMUL: compile-time enabled (runtime CPUID probe decides actual use)")
else()
  message(STATUS "PCLMUL: not available, using software carry-less multiply")
endif()

# Toeplitz 种子熵源：Windows 走 BCryptGenRandom，需显式链接 bcrypt
if(WIN32)
  target_link_libraries(qkdpp PUBLIC bcrypt)
endif()

enable_testing()
add_executable(qkdpp_tests tests/test_main.cpp)
target_link_libraries(qkdpp_tests PRIVATE qkdpp)
add_test(NAME unit COMMAND qkdpp_tests)

add_executable(qkdpp_bench bench/bench_toeplitz.cpp)
target_link_libraries(qkdpp_bench PRIVATE qkdpp)

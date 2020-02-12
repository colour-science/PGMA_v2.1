cmake_minimum_required(VERSION 3.15)
project(PGMA_v2_1 C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(PGMA_v2_1
        pgma_a_v2_1.c
        pgma_a_vgeom_v2_1.c
        pgma_gbcalc_v2_1.c
        pgma_gmas_v2_1.c
        pgma_test_v2_1.c
        pgma_v2_1.h)

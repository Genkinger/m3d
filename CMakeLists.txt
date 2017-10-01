cmake_minimum_required(VERSION 3.8)
project(m3d C)

set(CMAKE_C_STANDARD 99)

add_library(m3d SHARED m3d.c m3d.h)
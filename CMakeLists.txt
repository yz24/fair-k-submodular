cmake_minimum_required(VERSION 3.23)
project(ksubmodular)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(.)

add_library(IM
        mt19937ar.c
        jlog.h
        jlog.cpp
        tools.h
        algos.cpp
        picojson.h
        greedy.cpp
#        sensor.cpp
)
add_executable(ksubmodular mt19937ar.c
        jlog.h
        jlog.cpp
        tools.h
        picojson.h
        greedy.cpp
#        sensor.cpp
        )


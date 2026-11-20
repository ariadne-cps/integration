cmake_minimum_required(VERSION 3.15)
project(AriadneIntegrationTutorial LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(AriadneIntegration REQUIRED)

add_executable(tutorial_integration tutorial_integration.cpp)
target_link_libraries(tutorial_integration PRIVATE AriadneIntegration::ariadne-integration)

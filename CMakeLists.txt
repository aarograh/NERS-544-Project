cmake_minimum_required(VERSION 2.6)
project(NERS-544-Project)
add_subdirectory("${PROJECT_SOURCE_DIR}/src")
add_executable(NERS-544 ${SOURCE} ${HEADERS})

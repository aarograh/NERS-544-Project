# AUTHORS: Aaron Graham, Mike Jarrett
# PURPOSE: NERS 544 Course Project
# DATE   : April 3, 2015

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -ftrapv -Wextra -Wfloat-equal -Wundef -Wshadow -Wpointer-arith -Wcast-align -Wswitch-enum -Wconversion -Wunreachable-code -Wswitch-default"
)

set(CMake_CXX_FLAGS_RELEASE "-Ofast"
)

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG -g -Ofast"
)

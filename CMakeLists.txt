cmake_minimum_required(VERSION 3.20)
project(rsa_recon C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Compiler flags
add_compile_options(-Wall -Wextra -O2)

# Find required libraries
find_package(OpenSSL REQUIRED)
find_library(GMP_LIB gmp REQUIRED)

# Source files
set(SOURCES
        src/main.c
        src/common.c
        src/keygen.c
        src/decaySimulation.c
        src/initPhase.c
        src/heningerShacham.c
        src/heneckaMay.c
)

# Executable
add_executable(rsa_recon ${SOURCES})

# Include directories
target_include_directories(rsa_recon PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Link libraries
target_link_libraries(rsa_recon PRIVATE
        ${GMP_LIB}
        OpenSSL::SSL
        OpenSSL::Crypto
        m
)

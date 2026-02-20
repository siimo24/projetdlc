# ======================================================================
#  RSA Key Reconstruction — CMakeLists.txt
# ======================================================================

cmake_minimum_required(VERSION 3.20)
project(rsa_recon C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ---- Compiler flags ----
add_compile_options(-Wall -Wextra -O2)

# ---- Source files ----
set(SOURCES
        src/main.c
        src/common.c
        src/keygen.c
        src/decaySimulation.c
)

# ---- Executable ----
add_executable(rsa_recon ${SOURCES})

# ---- Include dirs ----
target_include_directories(rsa_recon PRIVATE include)

# ---- Dependencies ----
find_package(OpenSSL REQUIRED)

target_link_libraries(rsa_recon PRIVATE
        gmp
        OpenSSL::SSL
        OpenSSL::Crypto
        m
)

# ---- Output directory ----
set_target_properties(rsa_recon PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build
)

# ---- Custom targets (équivalents des convenience targets) ----
add_custom_target(run
        COMMAND rsa_recon --bits 1024 --preset cold --time 120 -v
        DEPENDS rsa_recon
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Running cold boot scenario (1024-bit)"
)

add_custom_target(test_quick
        COMMAND rsa_recon --bits 512 --preset room --time 5 -v
        DEPENDS rsa_recon
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Quick test with 512-bit key"
)

add_custom_target(simple
        COMMAND rsa_recon --bits 1024 --simple --known 0.30 --error 0.0
        DEPENDS rsa_recon
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Simple decay mode"
)

add_custom_target(frozen
        COMMAND rsa_recon --bits 1024 --preset frozen --time 3600 -v
        DEPENDS rsa_recon
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Frozen scenario (liquid nitrogen)"
)

add_custom_target(clean_all
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/keys
        COMMENT "Cleaning build artifacts"
)
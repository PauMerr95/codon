cmake_minimum_required(VERSION 3.16)
project(
	codon
	VERSION 0.1.0
	LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "CXX compiler ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Frontend variant: ${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}")
message(STATUS "Generator: ${CMAKE_GENERATOR}")

# -- This just makes the compilation output cleaner:
include(GNUInstallDirs)
set(OUTPUT_BASE ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_BASE}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${OUTPUT_BASE}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OUTPUT_BASE}/${CMAKE_INSTALL_LIBDIR})

# -- import compilation option (can be separately changed in cmake/compiler_options.cmake
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(compiler_options)

# -- specify libraries used --
add_library(codon_lib src/codon.cpp)


target_include_directories(codon_lib
	PUBLIC 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/plog/include
ncl)

# -- imported compiler options --
target_compile_options(codon_lib
	PRIVATE
	${codon_WARNING_FLAGS}
	${codon_BUILD_FLAGS}
)

# -- main executable / you can change the name here--
add_subdirectory(
    ${PROJECT_SOURCE_DIR}/include/plog
)

add_executable(codon src/main.cpp src/logging.cpp)

target_link_libraries(codon
	PRIVATE
    codon_lib
    plog::plog
)

target_compile_definitions(codon
    PRIVATE
    PLOG_FILE_NAME="${CMAKE_BINARY_DIR}/Log_codon.csv"
)

# -- Testing section --
enable_testing()

add_executable(codon_tests test/test_lib.cpp src/logging.cpp)

target_link_libraries(codon_tests
	PRIVATE
    codon_lib
    plog::plog
)

target_compile_definitions(codon_tests
    PRIVATE
    PLOG_FILE_NAME="${CMAKE_BINARY_DIR}/Log_codon_test.csv"
)

add_test(NAME codon_tests COMMAND codon_tests)

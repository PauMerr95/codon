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

message(STATUS "CXX compiler:       ${CMAKE_CXX_COMPILER}")
message(STATUS "CXX compiler ID:    ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Frontend variant:   ${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}")
message(STATUS "Generator:          ${CMAKE_GENERATOR}")

# -- This just makes the compilation output cleaner:
include(GNUInstallDirs)
set(OUTPUT_BASE ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_BASE}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${OUTPUT_BASE}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OUTPUT_BASE}/${CMAKE_INSTALL_LIBDIR})

# -- import compilation option (can be separately changed in cmake/compiler_options.cmake
# -- you can add your own or use the default flags
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(compiler_options OPTIONAL RESULT_VARIABLE compiler_opts_loaded_locally)

if(NOT compiler_opts_loaded_locally)
  message(WARNING "compiler_options.cmake not found. Reverting to default flags.")
  set(codon_WARNING_FLAGS "")
  set(codon_BUILD_FLAGS "")
endif()

# -- specify libraries used --
add_library(codon_lib src/codon.cpp src/seq.cpp)


target_include_directories(codon_lib
	PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/plog/include
)

# -- imported compiler options --
target_compile_options(codon_lib
	PRIVATE
	${codon_WARNING_FLAGS}
	${codon_BUILD_FLAGS}
)

# -- main executable / you can change the name here--
add_executable(codon src/main.cpp src/logging.cpp)

target_link_libraries(codon
	PRIVATE
    codon_lib
)

target_compile_definitions(codon
    PRIVATE
    PLOG_FILE_NAME="${CMAKE_BINARY_DIR}/Log_codon.csv"
)

# -- testing ---
option(ENABLE_COVERAGE "Enable coverage reporting" OFF)

if(ENABLE_COVERAGE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(STATUS "Coverage enabled")
        add_compile_options(--coverage -O0 -g)
        add_link_options(--coverage)
    else()
        message(WARNING "Coverage not supported for this compiler")
    endif()
endif()

include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.5.3
)

FetchContent_MakeAvailable(Catch2)

add_executable(testing
    test/test_instantiator.cpp
    test/test_codon.cpp
    test/test_seq.cpp
    src/logging.cpp)

target_link_libraries(testing
    PRIVATE
    codon_lib
    Catch2::Catch2WithMain)

include(CTest)
include(Catch)
catch_discover_tests(testing)

target_compile_definitions(testing
    PRIVATE
    PLOG_FILE_NAME="${CMAKE_BINARY_DIR}/Log_test_codon.csv"
)

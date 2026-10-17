cmake_minimum_required(VERSION 3.16)

project(fm_range_bwa_round1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_BENCHMARKS "Build standalone benchmark executables" OFF)

set(BWA_API_DIR "${CMAKE_CURRENT_SOURCE_DIR}/API" CACHE PATH "Directory containing bwt.h, bntseq.h, bwt.c, bntseq.c, utils.c")

if(NOT EXISTS "${BWA_API_DIR}/bwt.h")
    message(FATAL_ERROR "BWA API directory not found or missing bwt.h at ${BWA_API_DIR}")
endif()

# Build libbwa_bwt static library from C sources
add_library(bwa_bwt STATIC
    ${BWA_API_DIR}/bwt.c
    ${BWA_API_DIR}/bntseq.c
    ${BWA_API_DIR}/utils.c
)

target_include_directories(bwa_bwt PUBLIC ${BWA_API_DIR})

add_library(fm_range
    src/sa_range.cpp
    src/bwa_fmd_index.cpp
    src/mismatch_search.cpp
    src/exact_search.cpp
)

target_include_directories(fm_range
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${BWA_API_DIR}
)

target_link_libraries(fm_range PUBLIC bwa_bwt z)

add_executable(test_sa_range
    tests/test_sa_range.cpp
)

add_executable(test_bwa_extend_diagnostic
    tests/test_bwa_extend_diagnostic.cpp
)

add_executable(test_exact_search
    tests/test_exact_search.cpp
)

add_executable(test_bwa_fmd_exhaustive
    tests/test_bwa_fmd_exhaustive.cpp
)

add_executable(test_bwa_branch
    tests/test_bwa_branch.cpp
)

add_executable(test_bwa_locate
    tests/test_bwa_locate.cpp
)

add_executable(test_mismatch_search
    tests/test_mismatch_search.cpp
)

add_executable(test_fast_exact_search
    tests/test_fast_exact_search.cpp
)
target_link_libraries(test_sa_range PRIVATE fm_range)
target_link_libraries(test_bwa_extend_diagnostic PRIVATE fm_range)
target_link_libraries(test_exact_search PRIVATE fm_range)
target_link_libraries(test_bwa_fmd_exhaustive PRIVATE fm_range)
target_link_libraries(test_bwa_branch PRIVATE fm_range)
target_link_libraries(test_bwa_locate PRIVATE fm_range pthread m)
target_link_libraries(test_mismatch_search PRIVATE fm_range)
target_link_libraries(test_fast_exact_search PRIVATE fm_range)

if(BUILD_BENCHMARKS)
    add_executable(bench_exact_search
        benchmarks/bench_exact_search.cpp
    )
    target_link_libraries(bench_exact_search PRIVATE fm_range)
    target_compile_definitions(bench_exact_search PRIVATE
        BENCH_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
    )
endif()

enable_testing()

# The test programs need real BWA indexes (including sampled SA and, for
# locate(), reference metadata). Generate deterministic fixtures in the build
# tree rather than relying on generated files or absolute paths in the source
# checkout. The repository carries the BWA API sources, but not the `bwa`
# index-building executable.
find_program(BWA_INDEX_EXECUTABLE NAMES bwa)
if(NOT BWA_INDEX_EXECUTABLE)
    message(FATAL_ERROR
        "CTest fixtures require the BWA executable (bwa index) on PATH")
endif()

set(TEST_FIXTURE_ROOT "${CMAKE_CURRENT_BINARY_DIR}/test-fixtures")
set(EXACT_FIXTURE_DIR "${TEST_FIXTURE_ROOT}/exact")
set(LOCATE_FIXTURE_DIR "${TEST_FIXTURE_ROOT}/locate")
file(MAKE_DIRECTORY "${EXACT_FIXTURE_DIR}" "${LOCATE_FIXTURE_DIR}")

set(EXACT_FIXTURE_FASTA "${EXACT_FIXTURE_DIR}/reference.fa")
set(EXACT_FIXTURE_PREFIX "${EXACT_FIXTURE_DIR}/test")
file(WRITE "${EXACT_FIXTURE_FASTA}" ">test\nACGTTGCAACG\n")
execute_process(
    COMMAND "${BWA_INDEX_EXECUTABLE}" index
            -p "${EXACT_FIXTURE_PREFIX}" "${EXACT_FIXTURE_FASTA}"
    RESULT_VARIABLE EXACT_INDEX_RESULT
    OUTPUT_VARIABLE EXACT_INDEX_OUTPUT
    ERROR_VARIABLE EXACT_INDEX_ERROR
)
if(NOT EXACT_INDEX_RESULT EQUAL 0 OR
   NOT EXISTS "${EXACT_FIXTURE_PREFIX}.bwt" OR
   NOT EXISTS "${EXACT_FIXTURE_PREFIX}.sa")
    message(FATAL_ERROR
        "Failed to generate exact-search BWA test index with ${BWA_INDEX_EXECUTABLE}:\n${EXACT_INDEX_OUTPUT}${EXACT_INDEX_ERROR}")
endif()

set(LOCATE_FIXTURE_FASTA "${LOCATE_FIXTURE_DIR}/reference.fa")
set(LOCATE_FIXTURE_PREFIX "${LOCATE_FIXTURE_DIR}/test")
file(WRITE "${LOCATE_FIXTURE_FASTA}"
    ">ref1\nACGTTGCAACG\n>ref2\nTTACCGGTAA\n")
execute_process(
    COMMAND "${BWA_INDEX_EXECUTABLE}" index
            -p "${LOCATE_FIXTURE_PREFIX}" "${LOCATE_FIXTURE_FASTA}"
    RESULT_VARIABLE LOCATE_INDEX_RESULT
    OUTPUT_VARIABLE LOCATE_INDEX_OUTPUT
    ERROR_VARIABLE LOCATE_INDEX_ERROR
)
if(NOT LOCATE_INDEX_RESULT EQUAL 0 OR
   NOT EXISTS "${LOCATE_FIXTURE_PREFIX}.bwt" OR
   NOT EXISTS "${LOCATE_FIXTURE_PREFIX}.sa" OR
   NOT EXISTS "${LOCATE_FIXTURE_PREFIX}.ann" OR
   NOT EXISTS "${LOCATE_FIXTURE_PREFIX}.amb" OR
   NOT EXISTS "${LOCATE_FIXTURE_PREFIX}.pac")
    message(FATAL_ERROR
        "Failed to generate locate BWA test index with ${BWA_INDEX_EXECUTABLE}:\n${LOCATE_INDEX_OUTPUT}${LOCATE_INDEX_ERROR}")
endif()

add_test(NAME test_sa_range COMMAND test_sa_range)
add_test(NAME test_extend_diagnostics COMMAND test_bwa_extend_diagnostic "${EXACT_FIXTURE_PREFIX}.bwt")
add_test(NAME test_exact_search COMMAND test_exact_search "${EXACT_FIXTURE_PREFIX}.bwt" "${EXACT_FIXTURE_PREFIX}.sa")
add_test(NAME test_bwa_fmd_exhaustive COMMAND test_bwa_fmd_exhaustive "${EXACT_FIXTURE_PREFIX}" "${EXACT_FIXTURE_FASTA}")
add_test(NAME test_bwa_branch COMMAND test_bwa_branch "${EXACT_FIXTURE_PREFIX}.bwt" "${EXACT_FIXTURE_PREFIX}.sa")
add_test(NAME test_bwa_locate COMMAND test_bwa_locate "${LOCATE_FIXTURE_PREFIX}")
add_test(NAME test_mismatch_search COMMAND test_mismatch_search "${EXACT_FIXTURE_PREFIX}.bwt" "${EXACT_FIXTURE_PREFIX}.sa" "${EXACT_FIXTURE_FASTA}")
add_test(NAME test_fast_exact_search COMMAND test_fast_exact_search "${EXACT_FIXTURE_PREFIX}.bwt" "${EXACT_FIXTURE_PREFIX}.sa")

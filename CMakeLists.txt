cmake_minimum_required(VERSION 3.16)

project(fm_range_bwa_round1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
target_link_libraries(test_sa_range PRIVATE fm_range)
target_link_libraries(test_bwa_extend_diagnostic PRIVATE fm_range)
target_link_libraries(test_exact_search PRIVATE fm_range)
target_link_libraries(test_bwa_fmd_exhaustive PRIVATE fm_range)
target_link_libraries(test_bwa_branch PRIVATE fm_range)
target_link_libraries(test_bwa_locate PRIVATE fm_range pthread m)
target_link_libraries(test_mismatch_search PRIVATE fm_range)

enable_testing()
add_test(NAME test_sa_range COMMAND test_sa_range)
add_test(NAME test_extend_diagnostics COMMAND test_bwa_extend_diagnostic /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/minimal_bwt/index/test.fa.bwt)
add_test(NAME test_exact_search COMMAND test_exact_search /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/minimal_bwt/tests/tmp_exact_search/test.bwt /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/minimal_bwt/tests/tmp_exact_search/test.sa)
add_test(NAME test_bwa_fmd_exhaustive COMMAND test_bwa_fmd_exhaustive /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/minimal_bwt/tests/tmp_exact_search/test /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/minimal_bwt/tests/tmp_exact_search/reference.fa)
add_test(NAME test_bwa_branch COMMAND test_bwa_branch /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/minimal_bwt/tests/tmp_exact_search/test.bwt /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/minimal_bwt/tests/tmp_exact_search/test.sa)
add_test(NAME test_bwa_locate COMMAND test_bwa_locate /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/minimal_bwt/tests/tmp_locate/test)
add_test(NAME test_mismatch_search COMMAND test_mismatch_search /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/minimal_bwt/tests/tmp_exact_search/test.bwt /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/minimal_bwt/tests/tmp_exact_search/test.sa)

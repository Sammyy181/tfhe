cmake_minimum_required(VERSION 3.5)

# Build every test in this folder for every enabled FFT processor and
# register all of them with CTest.
#
#   cmake ../src -DENABLE_TESTS=ON -DCMAKE_BUILD_TYPE=debug
#   make -j
#   ctest --output-on-failure                 # everything except long-running tests
#   ctest -L unit                             # googletest unit tests only
#   ctest -L integration                      # standalone integration tests only
#   ctest -R multiplication                   # a single test by name
#   ctest -C long -L long                     # the 10M-gate long-run test

option(ENABLE_LONG_TESTS "Register test-long-run with CTest by default" OFF)

add_subdirectory(googletest EXCLUDE_FROM_ALL)

# The bundled googletest is old and trips newer compiler warnings, which the
# project-wide -Werror turns into errors. Don't let that block the tests.
foreach (GTEST_TARGET gtest gtest_main)
    if (TARGET ${GTEST_TARGET})
        target_compile_options(${GTEST_TARGET} PRIVATE -Wno-error)
    endif ()
endforeach ()

include_directories(${GTEST_INCLUDE_DIRS})

# googletest unit tests (all linked into one executable per FFT processor)
set(GOOGLETEST_SOURCES
        arithmetic_test.cpp
        lwe_test.cpp
        polynomial_test.cpp
        tlwe_test.cpp
        tgsw_test.cpp
        tlwe_fft_test.cpp
        tgsw_fft_test.cpp
        lwekeyswitch_test.cpp
        bootstrapping_test.cpp
        bootstrapping_test_fft.cpp
        io_test.cpp
        lagrangehalfc_test.cpp
        nussbaumer_test.cpp
        boots_gates_test.cpp
        fakes/lagrangehalfc.h
        fakes/lwe.h
        fakes/lwe-bootstrapping-fft.h
        fakes/lwe-bootstrapping.h
        fakes/tgsw-fft.h
        fakes/tgsw.h
        fakes/tlwe.h
        fakes/tlwe-fft.h
        fakes/lwe-keyswitch.h
        )

# Standalone integration tests: single source file, exit code 0 == pass.
set(CPP_ITESTS
        test-bootstrapping-fft
        test-bootstrapping-nussbaumer
        test-decomp-tgsw
        test-lwe
        test-multiplication
        test-tlwe
        test-gate-bootstrapping
        test-addition-boot
        )

set(C_ITESTS
        test-c-binding
        )

# Very slow (10M bootstrapped gates); always built, only run on request.
set(LONG_ITESTS
        test-long-run
        )

foreach (FFT_PROCESSOR IN LISTS FFT_PROCESSORS)

    if (FFT_PROCESSOR STREQUAL "fftw")
        set(RUNTIME_LIBS tfhe-fftw ${FFTW_LIBRARIES})
    else ()
        set(RUNTIME_LIBS tfhe-${FFT_PROCESSOR})
    endif ()

    # --- unit tests (googletest) ---
    set(UNIT_TARGET unittests-${FFT_PROCESSOR})
    add_executable(${UNIT_TARGET} ${GOOGLETEST_SOURCES} ${TFHE_HEADERS})
    target_link_libraries(${UNIT_TARGET} ${RUNTIME_LIBS} gtest gtest_main pthread)
    add_test(NAME ${UNIT_TARGET} COMMAND ${UNIT_TARGET})
    set_tests_properties(${UNIT_TARGET} PROPERTIES LABELS "unit;${FFT_PROCESSOR}")

    # --- C++ integration tests ---
    foreach (ITEST IN LISTS CPP_ITESTS LONG_ITESTS)
        set(ITEST_TARGET ${ITEST}-${FFT_PROCESSOR})
        add_executable(${ITEST_TARGET} ${ITEST}.cpp ${TFHE_HEADERS})
        target_link_libraries(${ITEST_TARGET} ${RUNTIME_LIBS})
    endforeach ()

    foreach (ITEST IN LISTS CPP_ITESTS)
        add_test(NAME ${ITEST}-${FFT_PROCESSOR} COMMAND ${ITEST}-${FFT_PROCESSOR})
        set_tests_properties(${ITEST}-${FFT_PROCESSOR} PROPERTIES
                LABELS "integration;${FFT_PROCESSOR}")
    endforeach ()

    # --- C integration tests ---
    foreach (ITEST IN LISTS C_ITESTS)
        set(ITEST_TARGET ${ITEST}-${FFT_PROCESSOR})
        add_executable(${ITEST_TARGET} ${ITEST}.c ${TFHE_HEADERS})
        target_link_libraries(${ITEST_TARGET} ${RUNTIME_LIBS})
        add_test(NAME ${ITEST_TARGET} COMMAND ${ITEST_TARGET})
        set_tests_properties(${ITEST_TARGET} PROPERTIES
                LABELS "integration;c;${FFT_PROCESSOR}")
    endforeach ()

    # --- long-running tests ---
    foreach (ITEST IN LISTS LONG_ITESTS)
        set(ITEST_TARGET ${ITEST}-${FFT_PROCESSOR})
        if (ENABLE_LONG_TESTS)
            add_test(NAME ${ITEST_TARGET} COMMAND ${ITEST_TARGET})
        else ()
            # only runs with `ctest -C long`
            add_test(NAME ${ITEST_TARGET} CONFIGURATIONS long COMMAND ${ITEST_TARGET})
        endif ()
        set_tests_properties(${ITEST_TARGET} PROPERTIES
                LABELS "long;${FFT_PROCESSOR}" TIMEOUT 0)
    endforeach ()

endforeach ()

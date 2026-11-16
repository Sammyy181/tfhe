# CMake generated Testfile for 
# Source directory: /Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test
# Build directory: /Users/samarthg/Desktop/Sem_7/BTP/tfhe/build-tests/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(unittests-nayuki-portable "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/build-tests/test/unittests-nayuki-portable")
set_tests_properties(unittests-nayuki-portable PROPERTIES  LABELS "unit;nayuki-portable" _BACKTRACE_TRIPLES "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;88;add_test;/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;0;")
add_test(test-bootstrapping-fft-nayuki-portable "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/build-tests/test/test-bootstrapping-fft-nayuki-portable")
set_tests_properties(test-bootstrapping-fft-nayuki-portable PROPERTIES  LABELS "integration;nayuki-portable" _BACKTRACE_TRIPLES "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;99;add_test;/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;0;")
add_test(test-bootstrapping-nussbaumer-nayuki-portable "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/build-tests/test/test-bootstrapping-nussbaumer-nayuki-portable")
set_tests_properties(test-bootstrapping-nussbaumer-nayuki-portable PROPERTIES  LABELS "integration;nayuki-portable" _BACKTRACE_TRIPLES "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;99;add_test;/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;0;")
add_test(test-decomp-tgsw-nayuki-portable "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/build-tests/test/test-decomp-tgsw-nayuki-portable")
set_tests_properties(test-decomp-tgsw-nayuki-portable PROPERTIES  LABELS "integration;nayuki-portable" _BACKTRACE_TRIPLES "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;99;add_test;/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;0;")
add_test(test-lwe-nayuki-portable "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/build-tests/test/test-lwe-nayuki-portable")
set_tests_properties(test-lwe-nayuki-portable PROPERTIES  LABELS "integration;nayuki-portable" _BACKTRACE_TRIPLES "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;99;add_test;/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;0;")
add_test(test-multiplication-nayuki-portable "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/build-tests/test/test-multiplication-nayuki-portable")
set_tests_properties(test-multiplication-nayuki-portable PROPERTIES  LABELS "integration;nayuki-portable" _BACKTRACE_TRIPLES "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;99;add_test;/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;0;")
add_test(test-tlwe-nayuki-portable "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/build-tests/test/test-tlwe-nayuki-portable")
set_tests_properties(test-tlwe-nayuki-portable PROPERTIES  LABELS "integration;nayuki-portable" _BACKTRACE_TRIPLES "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;99;add_test;/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;0;")
add_test(test-gate-bootstrapping-nayuki-portable "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/build-tests/test/test-gate-bootstrapping-nayuki-portable")
set_tests_properties(test-gate-bootstrapping-nayuki-portable PROPERTIES  LABELS "integration;nayuki-portable" _BACKTRACE_TRIPLES "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;99;add_test;/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;0;")
add_test(test-addition-boot-nayuki-portable "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/build-tests/test/test-addition-boot-nayuki-portable")
set_tests_properties(test-addition-boot-nayuki-portable PROPERTIES  LABELS "integration;nayuki-portable" _BACKTRACE_TRIPLES "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;99;add_test;/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;0;")
add_test(test-c-binding-nayuki-portable "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/build-tests/test/test-c-binding-nayuki-portable")
set_tests_properties(test-c-binding-nayuki-portable PROPERTIES  LABELS "integration;c;nayuki-portable" _BACKTRACE_TRIPLES "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;109;add_test;/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;0;")
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Ll][Oo][Nn][Gg])$")
  add_test(test-long-run-nayuki-portable "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/build-tests/test/test-long-run-nayuki-portable")
  set_tests_properties(test-long-run-nayuki-portable PROPERTIES  LABELS "long;nayuki-portable" TIMEOUT "0" _BACKTRACE_TRIPLES "/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;121;add_test;/Users/samarthg/Desktop/Sem_7/BTP/tfhe/src/test/CMakeLists.txt;0;")
endif()
subdirs("googletest")

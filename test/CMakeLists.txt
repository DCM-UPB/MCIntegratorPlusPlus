include_directories(common/)
link_libraries(mci)
if (USE_COVERAGE)
  link_libraries(--coverage)
endif()

add_executable(check main.cpp)
add_executable(ut1.exe ut1/main.cpp)
add_executable(ut2.exe ut2/main.cpp)
add_executable(ut3.exe ut3/main.cpp)

add_test(ut1 ut1.exe)
add_test(ut2 ut2.exe)
add_test(ut3 ut3.exe)

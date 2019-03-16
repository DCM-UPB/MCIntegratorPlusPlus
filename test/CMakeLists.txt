link_libraries(mci)

add_executable(check main.cpp)
#add_executable(ut1.exe ut1/main.cpp)
add_executable(ut2.exe ut2/main.cpp)
add_executable(ut3.exe ut3/main.cpp)
add_executable(ut4.exe ut4/main.cpp)

#add_test(ut1 ut1.exe)
add_test(ut2 ut2.exe)
add_test(ut3 ut3.exe)
add_test(ut4 ut4.exe)

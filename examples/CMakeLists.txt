link_libraries(mci)
add_executable(ex1.exe ex1/main.cpp)
add_executable(ex2.exe ex2/main.cpp)
target_link_libraries(ex2.exe ${MPI_CXX_LIBRARIES})

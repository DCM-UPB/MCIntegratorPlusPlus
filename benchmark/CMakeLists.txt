include_directories(common/)
link_libraries(mci)

add_executable(bench_integrate_nmc bench_integrate_nmc/main.cpp)

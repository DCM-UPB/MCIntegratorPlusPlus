file(GLOB SOURCES "*.cpp")
add_library(mci SHARED ${SOURCES})
add_library(mci_static STATIC ${SOURCES})

if (MPI_FOUND)
  target_link_libraries(mci ${MPI_CXX_LIBRARIES})
  target_link_libraries(mci_static ${MPI_CXX_LIBRARIES})
endif()

#First check if we can compile with XML2
try_compile(FLEUR_USE_XML ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_XML.f90
LINK_LIBRARIES ${FLEUR_LIBRARIES}
            )

if (NOT FLEUR_USE_XML)
      find_package(LibXml2)
      set(CMAKE_C_FLAGS "-I${LIBXML2_INCLUDE_DIR}")
      if (LIBXML2_LIBRARIES)
          set(TEST_LIBRARIES ${FLEUR_LIBRARIES} ${LIBXML2_LIBRARIES})
      endif()
try_compile(FLEUR_USE_XML ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_XML.f90
	    LINK_LIBRARIES ${TEST_LIBRARIES}
            )
       if (FLEUR_USE_XML)
              set(FLEUR_LIBRARIES ${TEST_LIBRARIES} )
       endif()
endif()

#Try to simply add -lxml2
if (NOT FLEUR_USE_XML)
      set(TEST_LIBRARIES ${FLEUR_LIBRARIES} -lxml2)

try_compile(FLEUR_USE_XML ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_XML.f90
            LINK_LIBRARIES ${TEST_LIBRARIES}
            )
       if (FLEUR_USE_XML)
              set(FLEUR_LIBRARIES -lxml2 ${FLEUR_LIBRARIES})
       endif()
endif()


message("XML Library found for linking:${FLEUR_USE_XML}")

if (FLEUR_USE_XML)
   try_compile(FLEUR_USE_XML ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_XML.c
   CMAKE_FLAGS "-DCMAKE_C_LINK_EXECUTABLE='echo no linking'" LINK_LIBRARIES "-lxml2")
   if (NOT FLEUR_USE_XML)
      find_package(LibXml2)
      set(CMAKE_C_FLAGS "-I${LIBXML2_INCLUDE_DIR}")
      try_compile(FLEUR_USE_XML ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/cmake/tests/test_XML.c
     )
   endif()
endif()

message("XML Library found for C:${FLEUR_USE_XML}")


if (FLEUR_USE_XML)
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_XML")
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_XML")
endif()

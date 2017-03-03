#Set the compiler names
set(CMAKE_Fortran_COMPILER mpif90)
set(CMAKE_C_COMPILER mpicc)
#Add include pathes
#set(CMAKE_C_FLAGS " -I$ENV{XML2_ROOT}/include")
if (DEFINED ENV{MAGMA_ROOT})
set(FLEUR_Fortran_FLAGS " -I$ENV{MAGMA_ROOT}/include")
set(FLEUR_LIBRARIES "-L$ENV{MAGMA_ROOT}/lib;-lmagma")
endif()
set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_AIX")
set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_AIX")

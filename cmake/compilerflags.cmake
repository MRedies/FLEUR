option(FLEUR_USE_OPENMP "whether to use OpenMP" ON)
#cmake file to set compiler flags for some of the known compilers
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   message("Intel Fortran detected")
   set(FLEUR_PRECISION_OPTION "-r8")
   set(FLEUR_OPENMP_FLAG "-qopenmp")
   if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "13.0.0.0")
      set(FLEUR_WARN_MESSAGE "You are using an old version of the Intel Fortran Compiler. Most likely FLEUR will not be build sucessfully. Consider to upgrade your compiler.")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl -assume byterecl")
      set(FLEUR_OPENMP_FLAG "-openmp")
      set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_OLDINTEL")
      set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_OLDINTEL")
   elseif (${CMAKE_Fortran_COMPILER_VERSION} VERSION_LESS "14.1.0.0")
      set(FLEUR_WARN_MESSAGE "You are using an old version of the Intel Fortran Compiler. The execution of the fleur_MPI might fail. Consider to upgrade your compiler.")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl -assume byterecl")
      set(FLEUR_OPENMP_FLAG "-openmp")
   elseif (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "17.0.0.0")
      set(FLEUR_WARN_MESSAGE "You are using an old version of the Intel Fortran Compiler. The execution of the fleur_MPI might fail. Consider to upgrade your compiler.")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl -assume byterecl")
   else()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl -qopenmp -assume byterecl -no-wrap-margin")
   endif()
   set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -xHost -O2 -g")
   #set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -xMIC-AVX512 -O2")
   if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "19.0.0.0")
       set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -C -traceback -O0 -g -check uninit -check pointers -DCPP_DEBUG -warn all")
   else()
       set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check arg_temp_created,assume,bounds,contiguous,format,output_conversion,pointers,stack,uninit -traceback -O0 -g -check uninit -check pointers -DCPP_DEBUG")
   endif()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
   message("PGI Fortran detected")
   set(FLEUR_PRECISION_OPTION "-Mr8;-Mr8intrinsics")
   set(FLEUR_OPENMP_FLAG "-mp")
   set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "") #fix problem in cmake
   #CPU
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -mp")
   #GPU
   #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mcuda=cuda9.0,cc60 -Mcudalib=cublas")
   #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mcuda:kepler+ -ta:tesla:cuda7.5 -DUSE_STREAMS -DNUM_STREAMS=${N_STREAMS} -Minfo=accel -acc")
   #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mcuda:cuda9.0,cc70 -DUSE_STREAMS -DNUM_STREAMS=${N_STREAMS} -Minfo=accel -acc")
   #set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fast -O3")
   set(CMAKE_Fortran_FLAGS_RELEASE "-O1 ") # to prevent cmake from putting -fast which causes problems with PGI18.4
   set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -C -traceback -O0 -g -Mchkstk -Mchkptr -Ktrap=fp -DCPP_DEBUG")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "XL")
   message("IBM/BG Fortran detected")
   set(FLEUR_PRECISION_OPTION "-qrealsize=8")
   set(FLEUR_OPENMP_FLAG "-qsmp=omp")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qnosave -qarch=qp -qtune=qp -qfixed -qsuppress=1520-022 -qessl")
   set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O4   -qsuppress=1500-036")
   set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}  -O0 -g -DCPP_DEBUG")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I/bgsys/local/libxml2/include/libxml2")
   set(FLEUR_DEFINITIONS ${FLEUR_DEFINITIONS} "CPP_AIX")
   set(FLEUR_MPI_DEFINITIONS ${FLEUR_MPI_DEFINITIONS} "CPP_AIX")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   message("gfortran detected")
   set(FLEUR_PRECISION_OPTION "-fdefault-real-8;-fdefault-double-8")
   set(FLEUR_OPENMP_FLAG "-fopenmp")
   if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "6.1.0")
      message(FATAL_ERROR "Only modern versions of gfortran >6.3 will be able to compile FLEUR\nYou need to specify a different compiler.\nSee the docs at www.flapw.de.\n")
   endif()
   #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -Wno-missing-include-dirs -DCPP_IRAPPROX")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -Wno-missing-include-dirs -fno-sign-zero")
   set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O2 -g")
   set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ffpe-trap=invalid,zero -fdump-core -Wall -Wextra -Wno-array-temporaries  -fbacktrace -fcheck=all  -finit-real=nan -O0 -g -DCPP_DEBUG")
else()
   message("Unknown compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
endif()

string(REPLACE ";" " " CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FLEUR_PRECISION_OPTION}")

if (FLEUR_USE_OPENMP)
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FLEUR_OPENMP_FLAG}")
endif()

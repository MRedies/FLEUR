hdf5_version=1.10.2
if [ ! -r CMake-hdf5-${hdf5_version} ]
then
    #Get the file with the code
    curl -LO "http://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-${hdf5_version}/src/CMake-hdf5-${hdf5_version}.tar.gz"
    tar xzf CMake-hdf5-${hdf5_version}.tar.gz
    cd CMake-hdf5-${hdf5_version}
    #copy options.cmake to adjust settings for compilation
    cp ../HDF5options.cmake .
    #Compile&test (This will take a while)
    ./build-unix.sh
    #Do make install 
    cd build
    make install
fi
#Store the installation location
FLEUR_LIBDIR="$PWD/HDF5/$hdf_version/lib $FLEUR_LIBDIR"
FLEUR_INCLUDEDIR="$PWD/HDF5/$hdf_version/include/static $FLEUR_INCLUDEDIR"

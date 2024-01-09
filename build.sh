rm -rf bin
mkdir build
cd build

cmake -D CMAKE_INSTALL_PREFIX=/Users/kakeru/Documents/Projects/FEM-ParallelFluidAnalisis/bin \
      ..

make && make install

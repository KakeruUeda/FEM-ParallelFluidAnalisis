mkdir build
cd build
cmake -D CMAKE_INSTALL_PREFIX=/Users/kakeru/Documents/Projects/cfdMPI/bin \
      ..

      # -D TP_DIR=/usr/local/TextParser \
      # -D enable_GLOG=ON \
      # -D CMAKE_INSTALL_PREFIX=/home/totani/bin/topologyOptim \
      # -D TP_DIR=/home/totani/lib/TextParser-1.8.5 \
      # -D GLPK_DIR=/home/totani/lib/glpk-4.64 \
      # -D EIGEN_DIR=/home/totani/lib/eigen-3.3.4 \
      # -D enable_GLOG=ON \
      # -D GLOG_DIR=/home/totani/lib/glog \

make && make install
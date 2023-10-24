cmake -P distclean.cmake
cmake -B build -S . -DCMAKE_BUILD_TYPE=Debug
cmake --build build

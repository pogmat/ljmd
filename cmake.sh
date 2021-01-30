#!/bin/bash

rm -r ./build/

mkdir build

cd build

echo "cmake .."
cmake -DENABLE_TESTING=ON ..

echo "cmake --build ."
cmake --build .

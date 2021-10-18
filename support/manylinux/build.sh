#!/bin/bash
cd /src
mkdir build
cd build
ls ..
cmake ..
make all -j2
ctest --no-compress-output -T test
make package


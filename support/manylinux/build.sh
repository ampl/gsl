#!/bin/bash
cd /src
mkdir build
chmod 777 build
cd build
cmake ..
make all -j2
ctest --no-compress-output -T test
make package


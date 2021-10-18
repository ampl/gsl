#!/bin/bash
cd `dirname $0`

IMG=mapgccv/manylinux1_x86_64:cmake-3.12.4
docker run -v `pwd`:/shared -v `pwd`/../../:/src --rm $IMG /shared/build.sh

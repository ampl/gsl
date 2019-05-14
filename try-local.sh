#!/bin/bash
# Usage example: ./try-local.sh -b gsl-trusty64 -b gsl-win64
buildbot try --connect=pb --master=127.0.0.1:5555 --username=gsl --passwd=gsl --vc=git $@



#!/bin/bash

TESTNAMES=`ls *.cpp`

for TEST in $TESTNAMES
do
    TEST=${TEST%%.*}
    echo "Processing $TEST ..."
    make TESTS=$TEST
done
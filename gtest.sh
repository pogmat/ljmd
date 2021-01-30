#!/bin/bash

if cd build/tests 2> /dev/null
then
    ctest
else
    echo "Please, compile first"
fi

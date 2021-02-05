#!/bin/bash

FILE_DIR=$1

cat << EOF
108
39.948
0.2379
3.405
8.5
17.1580
${FILE_DIR}/argon_108.rest
/dev/null
argon_108.dat
1000
5.0
100
EOF

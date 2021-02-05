#!/bin/bash

BINARY=$1
INPUT_DIR=$2

${INPUT_DIR}/input_files/argon_108.sh ${INPUT_DIR}/input_files | ${BINARY}

head -10 argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$1,$2,$3,$4);}'> a.dat
head -10 ${INPUT_DIR}/references/argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$1,$2,$3,$4);}'> b.dat
SUCCESS=0
cmp a.dat b.dat && SUCCESS=1

rm -f a.dat b.dat argon_108.dat

if [ ${SUCCESS} -ne 1 ]
then
    exit 1
fi

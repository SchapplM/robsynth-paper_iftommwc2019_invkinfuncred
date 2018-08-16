#!/bin/bash
# cd to this directory
for f in `find . -name *.txt`; do
  # sed -i "s/\\it//g" $f
  sed -i "s/S/S_/g" $f
  sed -i "s/C/C_/g" $f
done

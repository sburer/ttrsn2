#!/bin/bash

MATLAB="/Applications/MATLAB_R2022b.app/bin/matlab -nodesktop -nodisplay -nojvm"

cd mat
# for filename in $(ls instance_10_*.mat); do
for filename in $(ls instance_5_*.mat); do
  filename="${filename%.*}"
  echo "$filename"
  sed -e "s/INSERT_FILE_NAME_HERE/$filename.mat/g" ../run_test_template.m > $filename.m
  $MATLAB < "$filename".m > ../results/"$filename"_out.txt
  rm -rf $filename.m
done
cd ..


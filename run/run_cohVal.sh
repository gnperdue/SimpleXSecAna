#!/bin/sh

filelist="coh_valid.root coh_valid_15deg.csv coh_valid_30deg.csv coh_valid_allangles.csv"

# CC, 0.5, 10005
../bin/cohValidation -f gntp.10005.ghep.root
for file in $filelist
do
    mv ../hist/$file  ../hist/cc_0p5_$file
done

# NC, 0.5, 20005
../bin/cohValidation -f gntp.20005.ghep.root
for file in $filelist
do
    mv ../hist/$file  ../hist/nc_0p5_$file
done

# CC, 1.0, 10010
../bin/cohValidation -f gntp.10010.ghep.root
for file in $filelist
do
    mv ../hist/$file  ../hist/cc_1p0_$file
done

# NC, 1.0, 20010
../bin/cohValidation -f gntp.20010.ghep.root
for file in $filelist
do
    mv ../hist/$file  ../hist/nc_1p0_$file
done

# CC, 1.5, 10015
../bin/cohValidation -f gntp.10015.ghep.root
for file in $filelist
do
    mv ../hist/$file  ../hist/cc_1p5_$file
done

# NC, 1.5, 20015
../bin/cohValidation -f gntp.20015.ghep.root
for file in $filelist
do
    mv ../hist/$file  ../hist/nc_1p5_$file
done

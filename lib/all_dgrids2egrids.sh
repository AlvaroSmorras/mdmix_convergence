#!/bin/bash

dgrid_folder=$1
egrid_folder=$2

for dgrid in $(ls $dgrid_folder/*/*/*dx); do
    egrid=$(echo $dgrid | sed 's/'$dgrid_folder'/'$egrid_folder'/')
    egrid_subfolder=$(echo $egrid | cut -d'/' -f1,2,3)
    mkdir -p $egrid_subfolder
    python lib/dgrid2egrid.py $dgrid $egrid
done

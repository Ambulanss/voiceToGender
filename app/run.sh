#!/bin/sh
FILES=`ls resources/`
for file in $FILES 
do
    #echo $file
    python3 -W ignore app/inf132296_inf132206.py "resources/"$file
done

#!/bin/sh
FILES=`ls resources/`
for file in $FILES 
do
    echo $file
    python3 -W ignore app/main.py "resources/"$file
done

#!/bin/bash

#==============================================================
# This script return the columns of input file according to the
# key words.
#
# The first line of the input file indicate the column name.
#
# Inputs : $1 ---- INFILE.
#          $2 ---- Key words. Seperated by spaces.
#
# Output : Stdout
#
# Shule Yu
# Jun 20 2014
#==============================================================

NR=`echo "$2" | awk '{print NF}'`
NRL=`head -n 1 $1 | awk '{print NF}'`
for key in $2
do
    count=1

    for name in `head -n 1 $1`
    do
        if [ "${name}" = "${key}" ]
        then
            break
        fi

        count=$((count+1))
    done

    if [ "${count}" -le "${NRL}" ]
    then
        Index="${Index} ${count}"
    else
        echo "In Field Searching: No such field: ${key}" 1>&2
        echo "The first line in $1 is:" 1>&2
        head -n 1 $1 1>&2
        exit 1
    fi
done

awk -v A="${Index}" -v N=${NR} -v ORS=" " 'NR>1 {split(A,a);for (i=1;i<=N;i++) print $a[i];printf "\n"}' "$1"

exit 0

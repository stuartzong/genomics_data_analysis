#! /bin/bash
#
# This script is used to rename a patient in a pipeline results.
# need to be aggressively tested
#

while read line; do
    # get mapping keys
    dir=$(echo $line | awk '{print $1}'|tr -d "\n\r")
    old=$(echo $line | awk '{print $2}'|tr -d "\n\r")
    new=$(echo $line | awk '{print $3}'|tr -d "\n\r")
    echo "working on: $dir, $old, $new"
    cd $dir
    
    echo "replace file conten for: $dir"
    # find [^.]* -type f | xargs sed -i -s 's/TARGET-10-SJMPAL043775/xxxxx/g'
    find $dir -type f | xargs sed -i -s "s/$old/$new/g"
    
    echo "rename all files for: $dir"
    cd  $dir
    for i in $(find $dir -name "*$old*"); do echo "change file name for $i";mv -v $i $(echo $i | sed -s "s/$old/$new/"); done

    echo "change directory names"
    dir="/projects/trans_scratch/validations/genome-validator/NCI_ALL_MPAL/test"
    for i in $(find $dir -name "*$old*"); do mv -v $i $(echo $i | sed -s "s/$old/$new/"); done


done </projects/trans_scratch/validations/genome-validator/NCI_ALL_MPAL/test/c.tmp



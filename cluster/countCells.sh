#!/usr/bin/env bash

source /projects/sysbio/users/cellAtlas/env/bin/activate
out=/projects/sysbio/users/cellAtlas/data/dataset.meta/cellCounts.tsv
python_script=/projects/sysbio/users/cellAtlas/cluster/countCells.py

if [ -f $FILE ]; then 
    rm -f $FILE
fi

for file in `ls /projects/sysbio/users/cellAtlas/data/scanpyObj/*| grep .*.h5ad`; do
    python3 $python_script -i $file >> $out
done

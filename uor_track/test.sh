#!/bin/bash
lev=(850 500 250)
OUTDIR=/home/ys17-23/Extension2/renql/

cd ${OUTDIR} #match${suffix}

for file in ${lev[@]} ; do
    echo $file
done


#!/usr/bin/env bash

n=$1
n=${n%%.*}
n=${n##*/}
echo $n,`awk -F "," '{print $5}' $1 | awk -F "," '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'`,`awk -F "," '{print $5}' $1 | awk -F "," '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { print sqrt(mean2 / NR) / 10; }'`


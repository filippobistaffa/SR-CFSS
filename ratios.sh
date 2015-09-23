#!/usr/bin/env bash

sort -t "," -k 3,3 -n -o $1 $1
sort -t "," -k 3,3 -n -o $2 $2
n=$1
n=${n%%-*}
n=${n##*/}
#join -j 3 -t "," $1 $2 | awk -F "," '{print $5/$9}'
echo $n,`join -j 3 -t "," $1 $2 | awk -F "," '{print $5/$9}' | awk -F "," '{ sum += $1; n++ } END { if (n > 0) print sum / n; }'`,`join -j 3 -t "," $1 $2 | awk -F "," '{print $5/$9}' | awk -F "," '{delta = $1 - avg; avg += delta / NR; mean2 += delta * ($1 - avg); } END { print sqrt(mean2 / NR); }'`



#! /bin/bash

NODE_NAME="xeon_e5_2630_v3_$(uname -n | cut -d. -f1)_$(date +%s)"

for i in `seq 1 $1`;
do
	../../../orcc -v -e $2 2> ${NODE_NAME}.stderr
done

mkdir -p $NODE_NAME

./db2csv.py
mv results.* search_space.* *.stderr *.log $NODE_NAME

#!/bin/bash

codes="reference.x reference-vec.x increase-vec.x increase-vl.x flex-datatype.x"

echo -e "version\ttime_per_iteration"
for code in $codes
do
	OUT=`./$code`
	time=`echo $OUT | grep Micro | awk '{print $NF}'`
	echo -e "$code\t$time"
done


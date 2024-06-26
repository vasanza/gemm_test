#!/bin/bash

codes="test0.x test2.x test3.x test4.x test5.x version0.x reference-vec.x increase-vec.x increase-vl.x flex-datatype.x"

echo -e "version\ttime_per_iteration"
for code in $codes
do
	OUT=`./$code`
	time=`echo $OUT | grep Micro | awk '{print $NF}'`
	echo -e "$code\t$time"
done


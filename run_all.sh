#!/bin/bash

codes="test.x sgemm.x util.x ref.x"

echo -e "version\ttime_per_iteration"
for code in $codes
do
	OUT=`./$code`
	time=`echo $OUT | grep Micro | awk '{print $NF}'`
	echo -e "$code\t$time"
done

